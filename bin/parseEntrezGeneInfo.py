'''
Script to parse an Entrez GENE_INFO file to retrieve identifier, and gene names

gene_info file content description (from ftp://ftp.ncbi.nih.gov/gene/DATA/README):

gene_info                                       recalculated daily
---------------------------------------------------------------------------
           tab-delimited
           one line per GeneID
           Column header line is the first line in the file.
           Note: subsets of gene_info are available in the DATA/GENE_INFO
                 directory (described later)
---------------------------------------------------------------------------

[0] tax_id:
           the unique identifier provided by NCBI Taxonomy
           for the species or strain/isolate

[1] GeneID:
           the unique identifier for a gene
           ASN1:  geneid

[2] Symbol:
           the default symbol for the gene
           ASN1:  gene->locus

[3] LocusTag:
           the LocusTag value
           ASN1:  gene->locus-tag

[4] Synonyms:
           bar-delimited set of unofficial symbols for the gene

[5] dbXrefs:
           bar-delimited set of identifiers in other databases
           for this gene.  The unit of the set is database:value.

[6] chromosome:
           the chromosome on which this gene is placed.
           for mitochondrial genomes, the value 'MT' is used.

[7] map location:
           the map location for this gene

[8] description:
           a descriptive name for this gene

[9] type of gene:
           the type assigned to the gene according to the list of options
           provided in http://www.ncbi.nlm.nih.gov/IEB/ToolBox/CPP_DOC/lxr/source/src/objects/entrezgene/entrezgene.asn

[10] Symbol from nomenclature authority:
            when not '-', indicates that this symbol is from a
            a nomenclature authority

[11] Full name from nomenclature authority:
            when not '-', indicates that this full name is from a
            a nomenclature authority

[12] Nomenclature status:
            when not '-', indicates the status of the name from the 
            nomenclature authority (O for official, I for interim)

[13] Other designations:
            pipe-delimited set of some alternate descriptions that
            have been assigned to a GeneID
            '-' indicates none is being reported.

[14] Modification date:
            the last date a gene record was updated, in YYYYMMDD format
            

Note that in the below description we use the Entrez terminology. There are alternative terminologies:
Note that there are many variations for describing gene names. For exampe
- Systematic name = gene symbol
- Common name = gene symbol (used mostly in yeast)
- Alias = synonym

This script produces the following output files:

    * <KEGG ID>_gene_names.tab: a tab delimited text file where each line has the following format: 

  gene ID<tab>symbol<tab>synonym1|synonym2|...|synonymN<newline>

    * <KEGG ID>_synonym_to_id.tab: a tab delimited text file where each line has the following format : 

  synonym<tab>gene ID 1|gene ID 2|...|gene ID N<newline>

  Note that an alias can map to several IDs, and that the common name to ID mapping is always added.

    * <KEGG ID>_id_to_common.tab: a tab delimited text file where each line has the following format: 

  gene ID<tab>symbol<newline>

    * <KEGG ID>_gene_ids.txt: a text file with one gene ID per line.
    * <KEGG ID>.map: mapping of all symbols and synonyms to gene ID's:

  symbol or synonyms<tab>gene ID 1|gene ID 2|...|gene ID N<newline>
  
    * <KEG ID>_gene_list.txt: a tab delimited text file where each line has the following format:
    
  index<tab>gene ID<newline>, the first index is 1   

  
To handle gene name aliases the following rules are applied:
1. A gene symbol must map to exactly one gene ID. 
2. A gene synonym can not be a gene symbol name.
3. A gene synonym must map to exactly one gene ID.
4. The taxonomy ID must match the organisms ID

If a symbol maps to multiple gene ID's, only the first ID is included. If a synonym violates rule 2, 3, or 4 
it is discarded.

In addition it is possible to specify that:
4. There must be a specific database cross reference for a gene (the symbol/systematic name or synonym/alias 
 is discarded if the rule is violated).
 
orgSpecificDB = {
                 'xla': 'Xenbase:',
                 'sce': 'SGD:',
                 'mmu': 'MGI:',
                 'hsa': 'HGNC:',
                 'ath': 'TAIR:',
                 'rno': 'RGD:',
                 'ddi': 'dictyBase:'} 
'''

"""
Taxonomy ID's
"""
org2taxID = {
                 'xla': 8355,
                 'mmu': 10090,
                 'sce': 559292,
                 'hsa': 9606,
                 'ath': 3702,
                 'rno': 10116,
                 'ddi': 44689,
                 'dme': 7227,            
                 'cel': 6239} 

"""
Print usage description
"""
def usage():
    print """
Usage: python parseEntrezGeneInfo.py [options]\n

Required:
    -i FILENAME    Specify gene_info FILENAME for the organism.
    -o DIRECTORY   Specify output DIRECTORY.
    -k KEGG_ID     Specify the KEGG ID for the organism.
    
Optional:
    -a FILE1,FILE2 Specify a comma separated list of additional alias files with format: alias<tab>entrez ID<newline>-
    -d DATABASE    Specify a database crossreference to filter genes against (e.g. HGNC: for human)
    -t TYPE        Specify a gene type to filter genes against (e.g. protein-encoding) 
    -h             Display command line options.

"""

"""
Parse command line arguments. See the usage() output for the currently supported command
line arguments. Note also that most options are specified in the troilkatt.xml configuration
file.

@param args: command line arguments (sys.argv[1:])
"""
def parseArgs(argv):
    import getopt    
    
    try:
        opts, args = getopt.getopt(argv, "hi:o:k:a:d:t:")
    except getopt.GetoptError, err:
        # print help information and exit:
        print str(err)
        usage()
        sys.exit(2)
    
    args = {}
    for o, a in opts:
        if o == "-i":
            args['geneInfoFile'] = a
        elif o == "-o":
            args['outputDirectory'] = a
        elif o == "-k":
            args['keggID'] = a
        elif o == '-a':
            args['synonymFiles'] = a.split(',')
        elif o == '-d':
            args['dbXrefFilter'] = a
        elif o == '-t':
            args['type'] = a
        elif o in ("-h", "--help"):
            usage()
            sys.exit()        
        else:
            raise Exception("Unhandled option")
    
    if not args.has_key('geneInfoFile'):
        print 'Gene info file is not specified'
        usage()
        sys.exit(2)
    if not args.has_key('outputDirectory'):
        print 'Output directory is not specified'
        usage()
        sys.exit(2)
    if not args.has_key('keggID'):
        print 'KEGG ID is not specified'
        usage()
        sys.exit(2)
            
    return args

if __name__ == '__main__':
    import sys, os    
    
    args = parseArgs(sys.argv[1:])
    
    # Gene ID to symbol mapping
    id2symbol = {}
    # Gene symbol to ID mapping
    symbol2id = {}
    # In case of duplicate gene symbols, all the duplicate gene IDs in the alias files must be remapped
    duplicateIDRemap = {} # duplicate ID -> valid ID
    # Gene synonym to ID's mapping.
    synonym2ids = {}
    
    dbXrefFilter = None
    if args.has_key('dbXrefFilter'):
        dbXrefFilter = args['dbXrefFilter'].upper()
    
    if not org2taxID.has_key(args['keggID']):    
        raise Exception("Taxonomy ID not found for organism: %s" % (args['keggID']))
    orgTaxID = org2taxID[args['keggID']]
    
    typeFilter = None
    if args.has_key('type'):
        typeFilter = args['type'].upper()
    
    #
    # Read and parse gene_info file
    #
    inputFile = open(args['geneInfoFile'])
    geneNamesFile = open(os.path.join(args['outputDirectory'], args['keggID'] + '_gene_names.tab'), 'w')
    lineCnt = 0
    taxIDMismatch = 0
    duplicateSymbols = 0
    while 1:
        line = inputFile.readline()        
        if line == '':
            break
        
        lineCnt += 1
        # See above for file format
        parts = line.split('\t')
        
        # All are converted to upper case before mapping
        taxID = int(parts[0])
        geneID = parts[1].strip().upper()
        geneSymbol = parts[2].strip().upper()
        synonyms = []
        for s in parts[4].split('|'):
            if s != '-':
                synonyms.append(s.strip().upper())
        dbXrefs = parts[5].upper()
        dbSynonyms = []
        for d in dbXrefs.split('|'):
            dbSynonyms.append(d.strip())
        symbolSystematic = parts[10].strip().upper()
        systematic = parts[3].strip().upper()
        type = parts[9].strip().upper()        
        
        if taxID != orgTaxID:
            taxIDMismatch += 1
            continue
        
        if dbXrefFilter != None:        
            # Filter out genes without an official ID
            if dbXrefs.find(dbXrefFilter) == -1: # Not found
                continue
            
        if typeFilter != None:
            # Filter out genes that are not of the specified type
            if type.find(typeFilter) == -1: # Not found
                continue
                
        if id2symbol.has_key(geneID):
            raise Exception('Error: multiple entries with gene ID: %s' % (geneID))
                
        if symbol2id.has_key(geneSymbol):
            #raise Exception('Error: multiple entries with gene symbol: %s' % (geneSymbol))
            print 'Error on line %d: multiple entries with gene symbol: %s' % (lineCnt, geneSymbol)
            duplicateSymbols += 1
            duplicateIDRemap[geneID] = symbol2id[geneSymbol]
            continue
        
        geneNamesFile.write('%s\t%s\t%s\n' % (geneID, geneSymbol, parts[4]))
        id2symbol[geneID] = geneSymbol    
        symbol2id[geneSymbol] = geneID
        
        for s in synonyms:        
            if synonym2ids.has_key(s):
                synonym2ids[s].append(geneID)
            else:
                synonym2ids[s] = [geneID]

        for s in dbSynonyms:
            if synonym2ids.has_key(s):
                synonym2ids[s].append(geneID)
            else:
                synonym2ids[s] = [geneID]

        if symbolSystematic != '-' and symbolSystematic != '':
            if synonym2ids.has_key(symbolSystematic):
                synonym2ids[symbolSystematic].append(geneID)
            else:
                synonym2ids[symbolSystematic] = [geneID]

        if systematic != '-' and systematic != '':
            if synonym2ids.has_key(systematic):
                synonym2ids[systematic].append(geneID)
            else:
                synonym2ids[systematic] = [geneID]
        
    inputFile.close()
    geneNamesFile.close()
    
    #
    # Read and parse alias files 
    #
    if args.has_key('synonymFiles'):
        for sfn in args['synonymFiles']:
            sf = open(sfn)
            while 1:
                line = sf.readline()
                if line == '':
                    break
        
                # Each line has a gene name and the corresponding Entrez ID's
                parts = line.split('\t')                
                geneName = parts[0].strip().upper()
                for geneID in parts[1:]:                
                    # All gene names are treaded as synonyms
                    if synonym2ids.has_key(geneName):
                        synonym2ids[geneName].append(geneID.strip().upper())
                    else:
                        synonym2ids[geneName] = [geneID.strip().upper()]
            sf.close()
                    
    #
    # Create output files
    #    
    mapFile = open(os.path.join(args['outputDirectory'], args['keggID'] + '.map'), 'w')
    synonym2IDFile = open(os.path.join(args['outputDirectory'], args['keggID'] + '_synonym_to_id.tab'), 'w')
    
    multipleMappings = 0
    mapped = 0
    sameAsSymbol = 0
    invalidID = 0
    for s in symbol2id:
        mapFile.write('%s\t%s\n' % (s, symbol2id[s]))
    
    
    for a in synonym2ids:        
        # Ignore alias with same name as a symbol
        if symbol2id.has_key(a):
            sameAsSymbol += 1
            continue
        # Ignore alias that maps to multiple gene IDs
        if len(synonym2ids[a]) > 1:
            multipleMappings += 1
            continue
        
        synonymID = synonym2ids[a][0]
        # Ignore alias that does not map to a valid entrez ID
        if not id2symbol.has_key(synonymID):
            if duplicateIDRemap.has_key(synonymID):
                synonymID = duplicateIDRemap[synonymID]
            else:
                invalidID += 1
                print 'Invalid gene ID for synonym %s: %s' % (a, synonym2ids[a][0]) 
                continue
            
        mapFile.write('%s\t%s\n' % (a, synonym2ids[a][0]))
        synonym2IDFile.write('%s\t%s\n' % (a, synonym2ids[a][0]))
        mapped += 1
         
    print 'Total unique gene IDs and symbols: %d of %d' % (len(symbol2id), lineCnt)
    print '\tDuplicate symbols: %d' % (duplicateSymbols)
    print '\tTaxonomy ID mismatch for: %s genes' % (taxIDMismatch)
         
    print 'Mapped %d synonyms of %d'  % (mapped, len(synonym2ids))    
    print '\tSame name as a symbol: %d' % (sameAsSymbol)
    print '\tMaps to multiple gene IDs: %d' % (multipleMappings)
    print '\tInvalid gene ID: %d' % (invalidID) 
    mapFile.close()
    synonym2IDFile.close()
    
    geneIDs = symbol2id.keys()
    geneIDs.sort()
        
    geneIDsFile = open(os.path.join(args['outputDirectory'], args['keggID'] + '_gene_ids.txt'), 'w')            
    geneListFile = open(os.path.join(args['outputDirectory'], args['keggID'] + '_gene_list.txt'), 'w')
    
    for i in range(len(geneIDs)):
        geneIDsFile.write('%s\n' % (geneIDs[i]))
        geneListFile.write('%d\t%s\n' % (i + 1, geneIDs[i]))
        
    geneIDsFile.close()
    geneListFile.close()
    
