#!/usr/bin/env ruby
# Copyright (C) 2011 Peter Koppstein peak@princeton.edu
# License: Creative Commons Attribution-NonCommerical-ShareAlike 3.0 Unported (CC BY-NC-SA 3.0)
# See http://creativecommons.org/licenses/by-nc/3.0/
# Attribution shall include the copyright notice above.

# For help: $0 -h

# 0.0.4: handle equivalent GSEs properly
# 0.0.5: report GSEs for which |A-B|>2 (but see Revision 0.1.16)
# 0.0.6: "rule of 3" (excisions)
# 0.0.7: fix bug counting excisions; fix --gsm and --gse options
# 0.0.8: GSE = string
# 0.0.9: GSM dictionary
# 0.0.10: --excisions EXCISIONS ; handle missing entries from GSM dictionary explicitly
# 0.0.11: bugfix for excisions (introduce all_gses)
# 0.0.12: optimize EXCISIONS file; give counts for nontrivial overlaps; only run second round of excisions if necessary
# 0.0.13: co-ordinate documentation with analyzeGSMs
# 0.0.14: merge
# 0.1.16: @options.minimum

# Requirements: standard library

require 'optparse'
require 'ostruct'

# Using "set" might slow things down:
# In gene_gene_disease_count.rb, using set increases execution time by a factor of 3 (10 secs vs 30 secs)!
require 'set'

class Set
  # Return true iff the intersection is nonempty.
  # Without using cache: 19s; using cache: 28s; 
  def intersect? enum
    if enum === Set and length <= enum.length
      each { |o| return true if enum.include?(o)  }
    else
      enum.each { |o| return true if include?(o)  }
    end
    false
  end

  # return "min" if intersection has at least "min" elements, 
  # otherwise return | self & enum |
  def intersectionIndicator enum, min
    count=0
    each { |o| count += 1 if enum.include?(o)
      return min if count >= min
    }
    count
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
end # class Set

# Idiosyncratic methods:

def stringify_gsms that
    that.to_a.sort.each_with_index {|x,i| 
      if @options.dictionary
        print "GSM", x, ": ", (@gsm_dictionary[x] ? @gsm_dictionary[x] : ""), "\n"
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
end # class Set


# q.v. gene_gene_disease_count.rb

# Unless otherwise specified, for any key, options.key is nil
@options = OpenStruct.new

@options.version = "0.0.16"

@options.gse = 0     # column number minus 1
@options.gsm = 1     # column number minus 1
#@options.input_delimiter = "\t"
#@options.output_delimiter = "\t"

@options.input_delimiter = " "
@options.output_delimiter = " "

@options.minimum = 3

bn = File.basename($0)

optionparser = OptionParser.new do |opts|

  opts.banner = "Usage: #{bn} [options] FILE"

  opts.separator("")
  opts.separator("Options:")

  opts.on("", "--gse COLUMN", "Column number of gse (first column is numbered 1; default is #{@options.gse})") do |n|
    @options.gse = n.to_i - 1 # COLUMN is the column number, options.gsm is the column index
  end

  opts.on("", "--gsm COLUMN", "Column number of gsm (first column is numbered 1; default is #{@options.gsm})") do |n|
    @options.gsm = n.to_i - 1 # COLUMN is the column number, options.gsm is the column index
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

  opts.on("-m", "--minimum MINIMUM", "The minimum number of remaining samples there must be in order for an excision to be admissible") do |v|
    @options.minimum = v.to_i
  end

  opts.on("-t", "--terse", "Terse mode") do |n|
    @options.terse = true
  end

  opts.on("", "--excisions EXCISIONS", "Filename or pathname of file for excision directives (see below)") do |file|
    @options.excisions=file
  end

  opts.on("-v", "--[no-]verbose", "Run verbosely") do |v|
    @options.verbose = v
  end

  opts.on_tail("-V", "--version", "Show version and exit") do
    puts @options.version
    exit
  end

  explanation = <<-EOS
This script can be used to analyze the pattern of sample overlaps
amongst a set of GEO series files.  Optionally, a set of excision
directives can be written to a file, as explained below under
"Excision Directives".

In the following:

 . gse denotes an identifier of a family.soft file (usually this is a
   numeric GSE identifier, or a decimal identifier of the form x.y
   representing the series for platform y in GSE x)
 . gsm denotes the integer component of a GSM identifier
 . gsms(gse) is the set of integer gsm values associated with gse.

Input
-----
Unless otherwise specified by the options, input consists of
tab-delimited records of the form:

gse<tab>gsm

Normally, gse is an integer, n, representing GSEn, or a decimal of the
form x.y representing GPLy in GSEx.

Analysis
--------
This script analyzes the pattern of intersections amongst gsms(gse)
for a set of gse values, with a view to identifying gse files
that can be ignored, and to identifying a sequence of excision
directives that would be sufficient to eliminate non-trivial
overlaps:

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


Excision Directives
-------------------
An excision directive is a line with one of these forms:

gse ALL
gse [gsm ....]

That is, an excision directive is a gse value followed either by ALL
or by a possibly empty sequence of "gsm" values.

The first form means that the entire GEO series family.soft file is
redundant and can be removed from consideration.

The second form identifies a set of samples to be excised from the
corresponding gse file.  If no gsm values are given, nothing needs
to be excised.

If all the directives are followed, then the resultant family.soft
files will have few if any non-trivial overlaps.

The directives are derived by applying the following procedure
iteratively, where "minimum" defaults to 2:

    Consider two family.soft files, A and B, and let:
    a = gsms(A), the set of samples in A
    b = gsms(B), the set of samples in B.

    If |a & b| > 1 and |a| >= |b| and |a - b| >= minimum,
    then excise a&b from A.

The directives are then obtained by consolidating the excision operations.


GSM Dictionary
--------------
If the "--sample-title DICT" option is specified, where DICT is a
valid GSM dictionary, then sample titles will be printed alongside GSM
identifiers.

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

    #{bn} --sample-title gse_gsm.dict --excisions excisions.txt gse_gsm.txt

    #{bn} --sample-title data/download/yeast/gse_gsm.dict data/sce/gse_gsm.txt

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

all_gses = gses.dup  # freeze

@redundant = {}
@subsumed  = {}
@excised   = {}

############################################### START OUTPUT

print "Records for #{ngses} distinct GSE identifiers have been read.\n" unless @options.terse

if ngses == 0
  exit
end


@excisions = nil

if @options.excisions
  begin
    @excisions = File.open(@options.excisions, "w")
  rescue
    print "#{bn}: ERROR - unable to open file #{@options.excisions} but proceeding anyway\n"
  end
end

# gses:Array, gsms:Set
# Returns [overlap, nontrivial]
def overlaps(gses, gsms, verbose=false)
  overlap    = 0
  nontrivial = 0
  ngses = gses.length
  # Computing the intersection is what takes time
  gses.each_with_index { |gse1, i|
    (i+1 ... ngses).each { |j|
      gse2 = gses[j]
      ind = gsms[gse1].intersectionIndicator( gsms[gse2], 2 )
      if ind > 0
        overlap += 1
        nontrivial += 1 if ind > 1
        print "#{gse1} ^ #{gse2} = #{gsms[gse1].intersection( gsms[gse2] ).length}\n" if verbose
      end
    }
  }
  [overlap, nontrivial]
end

overlapCounts=overlaps(gses, gsms)

print "There are #{overlapCounts[0]} distinct cases of pairwise overlap.\n"
print "There are #{overlapCounts[1]} distinct cases of nontrivial pairwise overlap.\n"

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
      print "Non-trivial equivalence classes of GSEs:\n"
      printed = true
    end
    c.stringify
    print "\n"
  end
}

if @excisions
  equiv.each { |c| 
    c.each_with_index { |gse,i|
      if i > 0
        @excised[gse] = true
      end
    }
  }
end


# Select a representative from each equivalence class
# TODO: should we collect the greatest value of the GSE id?
equiv.collect! { |c| c.first }

gses = equiv.to_a.sort

print "There are #{gses.length} distinct GSE files.\n"

nduplicates = ngses - gses.length

print "Henceforth, unless otherwise specified, statistics ignore the #{nduplicates} GSEs that have been identified as duplicates.\n"

overlapCounts=overlaps(gses, gsms)

print "There are #{overlapCounts[0]} distinct cases of pairwise overlap (ignoring the duplicate GSEs).\n"
print "There are #{overlapCounts[1]} distinct cases of nontrivial pairwise overlap (ignoring the duplicate GSEs).\n"

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

print "\nIn the following, \"GSE X: S,T, ... \" means gsms(X) is a superset of each of gsms(S), ...\n"
print "and \": SUBSUMED\" indicates that gsms(X) is equal to the union of these sets.\n"

subsets.keys.each { |gse|
  # Compute gsms[gse] - UNION<x in subsets[gse]>  gsms[x]
  diff = gsms[gse].subtractOver(subsets[gse]) {|x| gsms[x]}
  # i.e. diff = gsms[gse].dup ; subsets[gse].each { |x| diff.subtract( gsms[x] ) }

  print "GSE #{gse}: "
  subsets[gse].stringify
  if diff.length == 0
    print " : SUBSUMED\n" 
    @subsumed[gse]  = 1
    @redundant[gse] = 1
    @excised[gse] = true
  else
    print "\n"
  end
}

print "\nThe number of GSEs that are 'SUBSUMED': #{@redundant.keys.length}\n"

# Recalculate ignoring the redundant values of gse
overlapping = {}
overlap = 0

# subsets[gse] is the set of GSEs which are subsets of gse
subsets = {}
gses.each_with_index { |gse1, i|
  gses.each_with_index { |gse2, j|
    next if i == j or @subsumed[gse1] or @subsumed[gse2]
    if gsms[gse1].subset?( gsms[gse2] )
      @redundant[gse1]=1
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
print "   # non-trivial union:          #{@subsumed.keys.length}\n"
print "   # proper inclusion:           #{@redundant.keys.length - @subsumed.keys.length}\n"

print "Total number of redundant GSEs (excluding duplicates):  #{@redundant.keys.length}\n"
print "\n"
print "# remaining distinct cases of overlap:                  #{overlap}\n"
print "# GSEs that overlap with another.                       #{overlapping.keys.length}\n"

###############################################
print "\nThe following connected components analyses are based on nontrivial intersection.\n"

print "Part 1: Global Analysis\n"
gses = Set[* gses]
cc = gses.divide { |g1,g2|
  # ignore trivial overlaps
  gsms[g1].intersectionIndicator( gsms[g2], 2 ) > 1    # i.e. (gsms[g1] & gsms[g2]).length > 1
}

connected_components = Set[]
cc.each { |c| connected_components.add(c) if c.length > 1 }
print "The number of non-trivial connected components (ignoring GSEs identified as duplicates) is #{connected_components.size}\n"
p connected_components


########################################################################
print "Part 2: Ignoring Redundant GSEs\n"
gses = gses - @redundant.keys
print "Examining #{gses.size} gses of interest...\n"

# return the number of remaining non-trivial overlaps
def excise connected_components, gses, gsms
  excisions = 0
  connected_components.each { |c|
    print "Connected component with #{c.length} GSEs:\n"
    c.each { |gse1|
      c.each { |gse2|
        next if (gse1 == gse2) or
                (gsms[gse1].length < gsms[gse2].length) or
                ! gsms[gse1].intersect?( gsms[gse2] )

        diff = gsms[gse1] - gsms[gse2]
        if diff.length >= @options.minimum  # default: there must be at least 3 samples left
          inter = (gsms[gse1] & gsms[gse2])
          if @excisions
            if @excised[gse1] 
              @excised[gse1] = (@excised[gse1] + inter)
            else
              @excised[gse1] = inter
            end
          end
          print "gsms[#{gse1}]: #{gsms[gse1].length} -> #{diff.length}\n"
          if @options.verbose and @options.dictionary
            printf "removing these GSMs from #{gse1}:\n"
            stringify_gsms inter
          end
          gsms[gse1] = diff
          excisions += 1
        end
      }
    }
  }

  print "Number of excisions: #{excisions}\n"
  
  print "\nPairwise overlaps after excision:\n"
  overlapping = overlaps(gses.to_a, gsms, true)
  print "   total overlaps:      #{overlapping[0]}\n"
  print "   nontrivial overlaps: #{overlapping[1]}\n"
  overlapping[1]
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
  print "If |A&B|>1 and |A| > |B| and A-B has at least #{@options.minimum} elements then use A-B and B\n"

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
      print "n-way intersection (n={#c.length}) = "; stringify_gsms all ; print "\n"
    end
  }

  print "The number of n-way intersections requiring review is #{overlap}\n"

  if @options.overlaps 
    print "\nThree-way overlaps:\n"
    connected_components.threeWayOverlaps(gsms)
  end

  #####

  print "\nExcising:\n"
  nti=excise connected_components, gses, gsms

  if nti > 0
    print "\nExcising - second round:\n"
    excise connected_components, gses, gsms
  end
}  # criterion

if @excisions
  # For the sake of neatness, proceed in three passes:
  all_gses.each { |gse|
    case @excised[gse]
    when true
      @excisions.print gse," ALL\n"
    end
  }
  all_gses.each { |gse|
    case @excised[gse]
    when Set
      @excisions.print gse
      @excised[gse].each {|gsm| @excisions.print " #{gsm}"}
      @excisions.print "\n"
    end
  }
  all_gses.each { |gse|
    case @excised[gse]
    when nil
      @excisions.print gse,"\n"
    end
  }
  @excisions.close
end
