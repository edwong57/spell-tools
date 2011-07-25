#!/usr/bin/env ruby
# Copyright (C) 2011 Peter Koppstein peak@princeton.edu
# License: Creative Commons Attribution-NonCommerical-ShareAlike 3.0 Unported (CC BY-NC-SA 3.0)
# See http://creativecommons.org/licenses/by-nc/3.0/
# Attribution shall include the copyright notice above.

# For help: $0 -h

# 0.0.4: handle equivalent GSEs properly
# 0.0.5: report GSEs for which |A-B|>2
# 0.0.6: "rule of 3" (excisions)
# 0.0.7: fix bug counting excisions; fix --gsm and --gse options
# 0.0.8: GSE = string
# 0.0.9: GSM dictionary

# Requirements: standard library

require 'optparse'
require 'ostruct'

# Using "set" might slow things down:
# In gene_gene_disease_count.rb, using set increases execution time by a factor of 3 (10 secs vs 30 secs)!
require 'set'

class Set
  # return true iff the intersection is nonempty
  # without using cache: 19s; using cache: 28s; 
  def intersect? enum
    if enum === Set and length <= enum.length
      each { |o| return true if enum.include?(o)  }
    else
      enum.each { |o| return true if include?(o)  }
    end
    false
  end

  # Usage: set.subtractOver(y) { |x| gsms[x] }
  def subtractOver y
    ans=dup
    y.each { |x| 
      ans.subtract yield x
    }
    ans
  end

  def stringify
    to_a.sort.each_with_index {|x,i| print "," if i>0; print x }
  end
end

# Idiosyncratic methods:

def stringify_gsms that
    that.to_a.sort.each_with_index {|x,i| 
      if @options.dictionary
        print "GSM", x, ": ", @gsm_dictionary[x], "\n"
      else
        print "," if i>0
        print x 
      end
  }
end

class Set 
  # self = connected_components
  def threeWayOverlaps(gsms)
    each { |c|
      print "\nConnected component with #{c.length} GSEs:\n"
      c.each { |gse1|
        c.each { |gse2|
          next if (gse1 <= gse2)
          int2 = gsms[gse1].intersection(gsms[gse2])
          next if int2.length == 0
          c.each { |gse3|
            next if (gse2 <= gse3)
            int3 = int2.intersection(gsms[gse3])
            next if int3.length == 0
            print "|#{gse1} ^ #{gse2} ^ #{gse3}| = #{int3.length}:\n"
            print "   "; stringify_gsms int3; print "\n"
          }
        }
      }
    }
  end
end


# q.v. gene_gene_disease_count.rb

# Unless otherwise specified, for any key, options.key is nil
@options = OpenStruct.new

@options.version = "0.0.9"

@options.gse = 0     # column number minus 1
@options.gsm = 1     # column number minus 1
#@options.input_delimiter = "\t"
#@options.output_delimiter = "\t"

@options.input_delimiter = " "
@options.output_delimiter = " "

optionparser = OptionParser.new do |opts|
  bn = File.basename($0)
  opts.banner = "Usage: #{bn} [options] FILE"

  opts.separator("")
  opts.separator("Options:")

  opts.on("", "--gse COLUMN", "Column number of gse (first column is numbered 1; default is #{@options.gse})") do |n|
    @options.gse = n.to_i - 1 # COLUMN is the column number, options.gsm is the column index
  end

  opts.on("", "--gsm COLUMN", "Column number of gsm (first column is numbered 1; default is #{@options.gsm})") do |n|
    @options.gsm = n.to_i - 1 # COLUMN is the column number, options.gsm is the column index
  end

  opts.on("-t", "--terse", "Terse mode") do |n|
    @options.terse = true
  end

  opts.on("", "--overlaps", "Provide information about three-way overlaps for GSEs that are not ignored") do |n|
    @options.overlaps = true
  end

  opts.on("", "--inputdelimiter CHAR", "Field delimiter for input, where CHAR may be tab or space for tab or space (default is tab)") do |n|
    # print "inputdelimiter=#{ n == "\t" ? "tab" : n == " " ? "space" : n}\n"
    n=" " if n == "space"
    n="\t" if n == "tab"
    @options.input_delimiter = n
  end

  opts.on("", "--outputdelimiter SEP", "Separator string for output, where SEP may be tab or space for tab or space (default is tab)") do |n|
    n=" " if n == "space"
    n="\t" if n == "tab"
    @options.output_delimiter = n
  end

  opts.on("", "--sample-title DICTIONARY", "Filename or pathname of file of or for GSM dictionary (see below)") do |dict|
    @options.dictionary=dict
  end

  opts.on("-v", "--[no-]verbose", "Run verbosely") do |v|
    @options.verbose = v
  end

  opts.on_tail("-V", "--version", "Show version and exit") do
    puts @options.version
    exit
  end

  explanation = <<-EOS
Let gse be a numeric GSE identifier, and let gsms(gse) be the set of
numeric GSM identifiers associated with gse.

This script analyzes the pattern of intersections amongst gsms(gse)
for a set of gse values:

1) identical: if gsms(gse1) == gsms(gse2) then one of them can be
   ignored;

2) complete redundancy: if one GSE is effectively just the union of
   some others that are not identical to it, then it can be ignored;

3) inclusion: if gse1 is not completely redundant and if there is
   another value, gse2, of GSE such that: gsms(gse1) is strictly
   contained in gsms(gse2), then gse1 can be ignored;

Amongst all the remaining GSE files, define the relation:
 "gse1 ~ gse2 if |gsms(gse1) & gsms(gse2)| > 1"

This relation defines a set of connected components.  This script
lists the connected components with more than one member.

Unless otherwise specified by the options, input consists of
tab-delimited records with integer values for gse and gsm in the first
two columns.

If the "--sample-title DICT" option is specified, where DICT is a
valid GSM dictionary (as defined in the next paragraph), then sample
titles will be printed alongside GSM identifiers.

A GSM dictionary is a file with lines of the form:

SERIES<tab>SAMPLE<tab>sample title

e.g.

GSE10018	GSM252981	Yeast aging 1 generation microarray 1

Lines beginging with "#" are ignored.

Notes:
FILE may be specified as - for stdin.

If a gsm id or gse id on a line is less than or equal to 0, the
line is ignored (silently unless the verbose option is given).

Long option names may be abbreviated to the shortest unique prefix.

Examples:
    #{bn} --input space -

    #{bn} --sample-title gse_gsm.dict gse_gsm.txt

    #{bn} --sample-title /Genomics/Users/peak/spell-tools/data/download/yeast/gse_gsm.dict /Genomics/Users/peak/spell-tools/data/sce/gse_gsm.txt

EOS

  opts.on("-h", "--help", "print this help") do
    puts opts
    print "\n", explanation
    exit 0
  end
end

# Modify ARGV
begin
  optionparser.parse!
rescue OptionParser::ParseError => e
  $stderr.print e, "\n"
  exit 1
end

if ARGV.length < 1
  puts optionparser.banner
  puts "For help type #{File.basename($0)} -h"
  exit 0
end

def verbose s
  print s + "\n" if @options.verbose
end


gsms = {}
ok=0

if @options.dictionary
  begin
    fdict = File.new(@options.dictionary, "r")
  rescue
    $stderr.print "file #{@options.dictionary} not found\n"
    exit 1
  end
end

@gsm_dictionary = {}
count=0

if fdict
  begin
    fdict.each_line { |line|
      line.chomp!
      next if line.start_with?("#")
      count += 1
      parts = line.split( "\t" )
      gsm=parts[0].sub("GSM","").to_i
      @gsm_dictionary[gsm] = parts[1]
    }
    verbose "#{count} items read from GSM dictionary"
  rescue
    $stderr.print "error reading GSM dictionary after reading #{count} line(s)\n"
    exit 1
  end
end


file = ARGV[0]
if file == "-"
  fd = STDIN
else
  begin
    verbose "Reading from #{ARGV[0]}"
    fd = File.new(file, "r")
  rescue
    $stderr.print "file #{file} not found\n"
    exit 1
  end
end

# IO.foreach(file)
count=0
fd.each_line do |line|
  count += 1
  line.chomp!
  parts = line.split( @options.input_delimiter )
  gse=parts[ @options.gse ]
  gsm=parts[ @options.gsm ].to_i

  if gsm > 0
    ok += 1
    if gsms[gse] 
      gsms[gse].add gsm
    else
      gsms[gse] = Set[gsm]
    end
  else
    $stderr.puts "skipping line #{count}: #{line}"
  end
end

fd.close unless file == "-"

verbose "#{count} lines read, #{ok} valid pairs found" 

sep = @options.output_delimiter

# This sort is solely to make the output tidy. 
# (Sorting fixnum is very fast.)
gses = gsms.keys.sort
ngses = gses.length

############################################### START OUTPUT

print "Records for #{ngses} distinct GSE identifiers have been read.\n" unless @options.terse

if ngses <= 1
  exit
end


# gses:Array, gsms:Set
def overlaps(gses, gsms, verbose=false)
  overlap = 0
  ngses = gses.length
  # Computing the intersection is what takes time
  gses.each_with_index { |gse1, i|
    (i+1 ... ngses).each { |j|
      gse2 = gses[j]
      if gsms[gse1].intersect?( gsms[gse2] )
        print "#{gse1} ^ #{gse2} = #{gsms[gse1].intersection( gsms[gse2] ).length}\n" if verbose
        overlap += 1
      end
    }
  }
  overlap
end

print "There are #{overlaps(gses, gsms)} distinct cases of pairwise overlap.\n"

# Select one representative from each equivalence class
# Warning: Set[].#first returns nil 
equiv = gses.to_set.divide { |g1,g2|
  (gsms[g1].length == gsms[g2].length) and
  (gsms[g1] == gsms[g2])
}

printed = nil
equiv.each { |c| 
  if c.length > 1 
    if not printed
      print "Equivalence classses of GSEs:\n"
      printed = true
    end
    c.stringify
    print "\n"
  end
}

equiv.collect! { |c| c.first }

gses = equiv.to_a.sort

print "There are #{gses.length} distinct GSE files.\n"

nduplicates = ngses - gses.length

print "Henceforth, unless otherwise specified, statistics ignore the #{nduplicates} GSEs that have been identified as duplicates.\n"

print "There are #{overlaps(gses, gsms)} distinct cases of pairwise overlap (ignoring the duplicate GSEs).\n"

redundant = {}
subsumed={}

# subsets[gse] is the set of GSEs which are subsets of gse
subsets = {}
gses.each_with_index { |gse1, i|
  gses.each_with_index { |gse2, j|
    next if i == j
    if gsms[gse1].subset?( gsms[gse2] )
      if subsets[gse2]
        subsets[gse2].add(gse1)
      else
        subsets[gse2] = Set[gse1]
      end
    end
  }
}


print "\nSUBSUMPTION\n"

subsets.keys.each { |gse|
  # Compute gsms[gse] - UNION<x in subsets[gse]>  gsms[x]
  diff = gsms[gse].subtractOver(subsets[gse]) {|x| gsms[x]}
  # i.e. diff = gsms[gse].dup ; subsets[gse].each { |x| diff.subtract( gsms[x] ) }

  print "GSE #{gse}: "
  subsets[gse].stringify
  if diff.length == 0
    print " : SUBSUMED\n" 
    subsumed[gse]  = 1
    redundant[gse] = 1
  else
    print "\n"
  end
}

print "\nThe number of GSEs that are 'SUBSUMED': #{redundant.keys.length}\n"

# Recalculate ignoring the redundant values of gse
overlapping = {}
overlap = 0

# subsets[gse] is the set of GSEs which are subsets of gse
subsets = {}
gses.each_with_index { |gse1, i|
  gses.each_with_index { |gse2, j|
    next if i == j or subsumed[gse1] or subsumed[gse2]
    if gsms[gse1].subset?( gsms[gse2] )
      redundant[gse1]=1
    elsif (i < j) and (gsms[gse1].intersect? gsms[gse2])
      overlapping[gse1] = 1
      overlapping[gse2] = 1
      overlap += 1
    end
  }
}

print "\nSummary:\n"

print "Number of GSEs including duplicates:                    #{ngses}\n"

print "   # duplicates:                 #{nduplicates}\n"
print "   # non-trivial union:          #{subsumed.keys.length}\n"
print "   # proper inclusion:           #{redundant.keys.length - subsumed.keys.length}\n"

print "Total number of redundant GSEs (excluding duplicates):  #{redundant.keys.length}\n"
print "\n"
print "# remaining distinct cases of overlap:                  #{overlap}\n"
print "# GSEs that overlap with another.                       #{overlapping.keys.length}\n"

###############################################
print "\nThe following connected components analyses are based on nontrivial intersection.\n"

print "Part 1: Global Analysis\n"
gses = Set[* gses]
cc = gses.divide { |g1,g2|
  (gsms[g1] & gsms[g2]).length > 1 # ignore trivial overlaps
  # gsms[g1].intersect?( gsms[g2] )   if gsms[g1]
  }

connected_components = Set[]
cc.each { |c| connected_components.add(c) if c.length > 1 }
print "The number of non-trivial connected components (ignoring GSEs identified as duplicates) is #{connected_components.size}\n"
p connected_components


########################################################################
print "Part 2: Ignoring Redundant GSEs\n"
gses = gses - redundant.keys
print "Examining #{gses.size} gses of interest...\n"

def excise connected_components, gses, gsms
  excisions = 0
  connected_components.each { |c|
    print "Connected component with #{c.length} GSEs:\n"
    c.each { |gse1|
      c.each { |gse2|
        next if (gse1 == gse2) or
        (gsms[gse1].length < gsms[gse2].length) or
        ! gsms[gse1].intersect?(gsms[gse2])
        diff = gsms[gse1] - gsms[gse2]
        if diff.length > 2  # there must be at least 3 samples left
          print "gsms[#{gse1}]: #{gsms[gse1].length} -> #{diff.length}\n"
          if @options.verbose and @options.dictionary
            printf "removing these GSMs from #{gse1}:\n"
            stringify_gsms (gsms[gse1] & gsms[gse2])
          end
            gsms[gse1] = diff
          excisions += 1
        end
      }
    }
  }

  print "Number of excisions: #{excisions}\n"
  
  print "\nOverlaps after excision:\n"
  overlapping = overlaps(gses.to_a, gsms, true)
  print "After excision, there are #{overlapping} distinct cases of pairwise overlap.\n"
end
  


[1].each { |criterion|

  print "criterion = #{criterion}...\n"

  cc = gses.divide { |g1,g2|
    (gsms[g1] & gsms[g2]).length > criterion # ignore trivial overlaps
    # gsms[g1].intersect?( gsms[g2] )   if gsms[g1]
  }

  connected_components = Set[]
  cc.each { |c| connected_components.add(c) if c.length > 1 }
  print "number of connected components is #{connected_components.size}\n"
  p connected_components
  
  print "The number of GSEs in these connected components is #{ connected_components.flatten.length }\n"

  #####
  print "If |A&B|>1 and |A| > |B| and A-B has at least three elements then use A-B and B"

  connected_components.each { |c|
    print "\n"
    c.each { |gse1|
      c.each { |gse2|
        next if (gse1 == gse2) or
          (gsms[gse1].length < gsms[gse2].length) or
          ! gsms[gse1].intersect?(gsms[gse2]) 
        diff = gsms[gse1] - gsms[gse2]
        if diff.length > 2
          print "| #{gse1} - #{gse2} | = #{diff.length}\n"
        end
      }
    }
  }

  #####
  # n-way intersection
  overlap = 0
  connected_components.each { |c|
    print "\nConnected component with #{c.length} GSEs:\n  "
    c.stringify
    print "\n"
    all = nil  # n-way intersection
    c.each { |gse1|
      all = all.nil?  ?  gsms[gse1] : all.intersection(gsms[gse1])
      c.each { |gse2|
        next if gse1 >= gse2
        inter = gsms[gse1].intersection(gsms[gse2])
        if inter.length > 0
          print "#{gse1} ^ #{gse2} = "
               print "#{gsms[gse1].length} ^ #{gsms[gse2].length} = ", inter.length
               print "  :: "
          print "\n" if @options.dictionary
          stringify_gsms inter
          print "\n"
        end
      }
    }
    if all === Set and all.length > 0
      overlap += 1 if inter.length > 1
      print "all="; stringify_gsms all ; print "\n"
    end
  }

  print "The number of n-way intersections requiring review is #{overlap}\n"

  if @options.overlaps 
    print "\nThree-way overlaps:\n"
    connected_components.threeWayOverlaps(gsms)
  end

  #####

  print "\nExcising:\n"  
  excise connected_components, gses, gsms

    print "\nExcising - second round:\n"
  excise connected_components, gses, gsms
}  # criterion

