#!/usr/bin/ruby

require 'set'

def cancelout()
  $stderr.print("Received negative number; canceling...\n")
  Process.exit
  return
end

if ARGV.length < 1
  puts "ARGS:"
  puts "0 - GSE***_family.soft file"
  exit 0
end

inPlatTable = false
inSampTable = false
inSample = false

platIdIdxs = Hash.new			#Per-platform, keep which column has unique IDs
platNameIdxs = Hash.new		#Per-platform, keep which column has gene names

sampIdIdxs = Hash.new 		#Per-platform, which column of the sample table has unique IDs
sampValueIdxs = Hash.new 	#Per-platform, which column of the sample table has desired value

sampPlatform = Hash.new		#per-sample, which platform to use
sampTitle = Hash.new			#per-sample, title string for sample

platLines = nil
sampLines = nil
platID = nil
sampID = nil

missingValueCount = 0
zeroValueCount = 0

#Adding config script hack

if ARGV.length ==2
  $stderr.print("Will use parameters from config file\n")
  useConfig=true
  file = File.new(ARGV[1], "r")
  line=file.gets.to_i
  platIdIdxAll=line

  line=file.gets.to_i
  platNameIdxAll=line

  line=file.gets.to_i
  sampIdIdxAll=line

  line=file.gets.to_i
  sampValueIdxAll=line

  line=file.gets.to_i
  convertZerosToMV=line

  $stderr.printf("Platform ProbeId is in column %d\n", platIdIdxAll)
  $stderr.printf("Platform GeneName is in column %d\n", platNameIdxAll)
  $stderr.printf("Sample ProbeId is in column %d\n", sampIdIdxAll)
  $stderr.printf("Sample Value is in column %d\n", sampValueIdxAll)
  if convertZerosToMV==1
    $stderr.print("Will convert zeros to missing values\n")
  else
    $stderr.print("Will NOT convert zeros to missing values\n")
  end
else
  useConfig=false
end

#Do an initial pass through to determine information about platforms, samples, etc.
IO.foreach(ARGV[0]) do |line|
  line.chomp!

  if line.include?("^PLATFORM")
    platID = line.slice(12,line.length)
    platLines = Array.new
  end
  if line.include?("!platform_table_begin")
    inPlatTable = true
  elsif inPlatTable == true
    if line.include?("!platform_table_end")
      inPlatTable = false
    elsif platLines.length < 5
      platLines.push(line)
    end
    if platLines.length == 5 && platIdIdxs[platID] == nil
      if useConfig
        platIdIdxs[platID] = platIdIdxAll
        platNameIdxs[platID] = platNameIdxAll
      else
        $stderr.puts "Columns for platform " + platID
        numCols = platLines[0].split("\t").length
        for j in 0...numCols
          $stderr.print j.to_s
          for i in 0...platLines.length
            parts = platLines[i].split("\t")
            if parts.length > j
              $stderr.print "\t" + parts[j]
            else
              $stderr.print "\t--"
            end
          end
          $stderr.puts
        end

        $stderr.print "Which row contains the probe ID? "
        platIdIdxs[platID] = $stdin.gets.to_i
        # useful if called from batch file -- pjbradle
        # ADD END, USEFUL MSG
        if platIdIdxs[platID] < 0: cancelout() end
        $stderr.print "Which row contains the gene names? "
        platNameIdxs[platID] = $stdin.gets.to_i
        if platNameIdxs[platID] < 0: cancelout() end
      end
    end
  end

  if line.include?("^SAMPLE")
    inSample = true
    sampID = line.slice(10,line.length)
  elsif inSample == true
    if line.include?("!Sample_title")
      sampTitle[sampID] = line.slice(16,line.length)
    elsif line.include?("!Sample_platform_id")
      sampPlatform[sampID] = line.slice(22,line.length)
    elsif line.include?("!sample_table_begin")
      inSampTable = true
      sampLines = Array.new
    elsif inSampTable == true
      if line.include?("!sample_table_end")
        inSampTable = false
        inSample = false
      elsif sampLines.length < 5
        sampLines.push(line)
      end
      if sampLines.length == 5 && sampIdIdxs[sampPlatform[sampID]] == nil
        if useConfig
          sampIdIdxs[sampPlatform[sampID]] = sampIdIdxAll
          sampValueIdxs[sampPlatform[sampID]] = sampValueIdxAll
        else
          $stderr.puts "Columns for sample " + sampID + " using platform " + sampPlatform[sampID].to_s
          numCols = sampLines[0].split("\t").length
          for j in 0...numCols
            $stderr.print j.to_s
            for i in 0...sampLines.length
              parts = sampLines[i].split("\t")
              if parts.length > j
                $stderr.print "\t" + parts[j]
              else
                $stderr.print "\t--"
              end
            end
            $stderr.puts
          end

          $stderr.print "Which row contains the probe ID? "
          sampIdIdxs[sampPlatform[sampID]] = $stdin.gets.to_i
          if sampIdIdxs[sampPlatform[sampID]] < 0: cancelout() end
          $stderr.print "Which row contains the expression values? "
          sampValueIdxs[sampPlatform[sampID]] = $stdin.gets.to_i
          if sampValueIdxs[sampPlatform[sampID]] < 0: cancelout() end
        end
      end
      if sampIdIdxs[sampPlatform[sampID]] != nil
        parts = line.split("\t")
        val = parts[sampValueIdxs[sampPlatform[sampID]]]
        if val == ""
          missingValueCount += 1
        elsif val.to_f == 0
          zeroValueCount += 1
        end
      end
    end
  end
end
if useConfig
  zerosAsMissingVals=convertZerosToMV
else
  $stderr.puts "\nThere are " + missingValueCount.to_s + " missing values, and " + zeroValueCount.to_s + " exactly zero values"
  $stderr.print "Treat exact zeros as missing values (0=no, 1=yes)? "
  zerosAsMissingVals = $stdin.gets.to_i
end
$stderr.puts
$stderr.puts "Parsing file into pcl format..."

samples = Hash.new	#Hash from sampID to a Hash from geneName to value
platMaps = Hash.new	#Hash from platID to a Hash from uniqueID to geneName
idToName = nil
plat = nil
genes = Set.new
first = false

#On a second pass through, actually parse out everything
IO.foreach(ARGV[0]) do |line|
  line.chomp!

  if line.include?("^PLATFORM")
    platID = line.slice(12,line.length)
    platMaps[platID] = Hash.new
  end
  if line.include?("!platform_table_begin")
    inPlatTable = true
  elsif inPlatTable == true
    if line.include?("!platform_table_end")
      inPlatTable = false
    else
      parts = line.split("\t")
      platMaps[platID][parts[platIdIdxs[platID]]] = parts[platNameIdxs[platID]]
    end
  end

  if line.include?("^SAMPLE")
    inSample = true
    sampID = line.slice(10,line.length)
    samples[sampID] = Hash.new
    plat = sampPlatform[sampID]
    idToName = platMaps[plat]
  elsif inSample == true
    if line.include?("!sample_table_begin")
      inSampTable = true
      first = true
    elsif inSampTable == true
      if first == true
        first = false
      elsif line.include?("!sample_table_end")
        inSampTable = false
        inSample = false
        plat = nil
        idToName = nil
      else
        parts = line.split("\t")
        gname = idToName[parts[sampIdIdxs[plat]]]
        if gname != nil && gname != ""
          samples[sampID][gname] = parts[sampValueIdxs[plat]]
          genes.add(gname)
        else
          #$stderr.puts "WARNING: For sample " + sampID + ", " + parts[sampIdIdxs[plat]] + " was not present in the " + sampPlatform[sampID] + " platform table"
        end
      end
    end
  end
end

#Output the pcl file
sampArr = samples.keys.to_a.sort
print "YORF\tNAME\tGWEIGHT"
for i in 0...sampArr.length
  print "\t" + sampTitle[sampArr[i]]
end
puts
print "EWEIGHT\t\t"
for i in 0...sampArr.length
  print "\t1"
end
puts

for gene in genes.to_a.sort
  print gene + "\t" + gene + "\t1"
  for i in 0...sampArr.length
    val = samples[sampArr[i]][gene]
    if val == nil
      print "\t"
    else
      if zerosAsMissingVals == 1 && val.to_f == 0
        print "\t"
      else
        print "\t" + val
      end
    end
  end
  puts
end
