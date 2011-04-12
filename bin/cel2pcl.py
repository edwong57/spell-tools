import os
from troilkatt_script import TroilkattScript

"""
Convert a set of CEL files to a PCL file
"""
class Cel2Pcl(TroilkattScript):
    """
    arguments: 
    [0] Directory with R script and R libriaries
    [1] Organism code (hs = human, rn = rat, and so on)
        
    @param: see description for super-class
    """
    def __init__(self):
        TroilkattScript.__init__(self)
                
        argsParts = self.args.split(" ")
        if (len(argsParts) != 2):
            raise Exception('Invalid arguments: %s' % (self.args))
        self.rDir = argsParts[0]
        self.orgCode = argsParts[1]

    """
    Cleanup input directory by removing all files.

    @param inputFiles: list of files that should not be removed
    @return none
    """
    def cleanupInputDir(self, inputFiles=None):
        os.chdir(self.inputDir)
        dirContent = os.listdir('.')

        inputFileBasenames = []
        for fn in inputFiles:
            inputFileBasenames.append(os.path.basename(fn))
        
        for fn in dirContent:
            if fn in inputFileBasenames:
                continue
            else:
                os.remove(fn)
    
    """
    Script specific run function. A sub-class should implement this function.
    
    @param inputFiles: list of absolute filenames in the input directory
    @param inputObjects: dictionary with objects indexed by a script specific key (can be None).
    
    @return: dictionary with objects that should be saved in Hbase using the given keys (can be None).
    """
    def run(self, inputFiles, inputObjects):
        oldDir = os.getcwd()
        self.cleanupInputDir(inputFiles)
        #
        # 1. Unpack tar file
        #
        tarName = None
        for fn in inputFiles:
            if self.endsWith(fn, '.tar'):
                tarName = os.path.basename(fn).split(".")[0]
                os.chdir(self.inputDir)
                cmd = 'tar xvf %s > %s 2> %s' % (fn, 
                                                 os.path.join(self.logDir, os.path.basename(fn) + '.untar.output'),
                                                 os.path.join(self.logDir, os.path.basename(fn) + '.untar.error'))
                #print 'Execute: %s' % (cmd)
                if os.system(cmd) != 0:
                    print 'Unpack failed: %s' % (cmd)
                    self.logger.warning('Unpack failed: %s' % (cmd))
                    self.cleanupInputDir()
                    os.chdir(oldDir)
                    return None
                
        if tarName == None:
            raise Exception('RAW.tar file not found in input directory: ' % (self.inputDir))
        
        #
        # 2. Decompress CEL files
        #
        unpackedFiles = os.listdir(self.inputDir)
        for fn in unpackedFiles:
            if self.endsWith(fn.lower(), 'cel.gz'):
                cmd = 'gunzip -f %s > %s 2> %s' % (fn, 
                                                os.path.join(self.logDir, fn + '.gunzip.output'),
                                                os.path.join(self.logDir, fn + '.gunzip.error'))
                #print 'Execute: %s' % (cmd)
                if os.system(cmd) != 0:
                    print 'gunzip failed: %s' % (cmd)
                    self.logger.warning('gunzip failed: %s' % (cmd))
                    self.cleanupInputDir()
                    os.chdir(oldDir)
                    return None
                    
        
        #
        # 3. Run R script to create PCL file in input directory 
        #
        outputPrefix = tarName + '.tmp'
        rOutputPrefix = os.path.join(self.inputDir, outputPrefix)
        print 'Chdir to: %s' % (self.rDir)
        os.chdir(self.rDir)
        cmd = '/nhome/larsab/troilkatt/apps/bin/troilkatt-container 12 -1 ./R --no-save --args %s %s %s < ProcessCEL.R > %s 2> %s' % (self.inputDir, 
                                                                                                                                     self.orgCode,
                                                                                                                                     rOutputPrefix,
                                                                                                                                     os.path.join(self.logDir, 'R.output'),
                                                                                                                                     os.path.join(self.logDir, 'R.error'))
        
        
        print 'Execute: %s' % (cmd)
        if os.system(cmd) != 0:
            print 'R script failed'
            self.logger.warning('R script failed')
            self.cleanupInputDir()            
            os.chdir(oldDir)
            return None

        #
        # 4. Merge and convert partial files
        #      
        partFiles = os.listdir(self.inputDir)
        for f in partFiles:
            if f.find(outputPrefix) != -1:
                if f.find(".single") != -1 or f.find(".zero") != -1:
                    continue
                
                if f.find('platform') != -1:
                    newName = f.replace('.tmp.platform.', '_') + '.pcl'                
                else:
                    newName = f.replace('.tmp', '.pcl')

                fin = open(os.path.join(self.inputDir, f))
                fout = open(os.path.join(self.outputDir, newName), 'w')
                headerLine = fin.readline()
                # Add extra columns for additional gene ID and GWEIGHT
                fout.write('ENTREZ_ID\tENTREZ_ID\tGWEIGHT')
                headerParts = headerLine.split('\t')
                for h in headerParts[1:]: # Ignore empty first column
                    fout.write('\t' + h)
                # last column contains newline

                # Add EWEIGHT column
                fout.write('EWEIGHT\t\t1')
                for h in headerParts[1:]:
                    fout.write('\t1')
                fout.write('\n')

                # Fix gene IDs and add the two additional columns
                while 1:
                    l = fin.readline()
                    if l == '':
                        break

                    cols = l.split('\t')
                    if len(cols) < 2:
                        continue
                    fout.write('%s\t%s\t1' % (cols[0], cols[0]))
                    for c in cols:
                        fout.write('\t' + c.strip())
                    # last column has newline
                fin.close()
                fout.close()
                    
        self.cleanupInputDir()
        os.chdir(oldDir)
                
        return None
        
"""
Run a troilkatt script

Command line arguments: %s input-dir output-dir log-dir args, where
   input-dir      Directory containing files to process (or where to store downloaded files)
   output-dir     Directory where output files are stored.
   log-dir        Directory where logfiles are stored.
   args[0]        Directory with R script
   args[1]        Organism code used by R script

The decription in usage() has additional details. 
"""
if __name__ == '__main__':
    import os
    os.chdir('/nhome/larsab/skarntyde/troilkatt-java')
    s = Cel2Pcl()
    s.mainRun()
