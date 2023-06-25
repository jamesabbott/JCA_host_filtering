import pandas as pd
import pickle
import os.path
import tarfile
from ftplib import FTP

def download_file(data_dir, ftphost, ftpdir, ftpfile):

    """
    Download file via ftp if newer than existing file, storing remote timestamp
    Download file via ftp if newer than existing file, storing remote timestamp
    in filename.timestamp file

    arguments:
    none

    returns:
    none
    """

    ftp = FTP(ftphost)
    ftp.login()
    ftp.cwd(ftpdir)
    timestamp = ftp.sendcmd('MDTM '+ ftpfile)
    # Successful command return will be prefixed with '213' return code
    if timestamp.startswith('213 '):
        timestamp=timestamp.replace('213 ', '')
    else:
        print("Error retrieving timestamp of %s" % ftpfile)
        exit(1)

    ts_file = ('%s/%s.timestamp' % (data_dir,ftpfile))
    if os.path.isfile(ts_file):
        with open(ts_file, 'r') as ts:
            old_ts = ts.read()
    else:
        old_ts = ''

    if timestamp != old_ts:
        print("Downloading new version of %s..." % ftpfile )
        outfile = open('%s/%s' % (data_dir,ftpfile), 'wb')

        ftp.retrbinary('RETR %s' % ftpfile, outfile.write)
        outfile.close()

        # save timestamp of downloaded file for future comparions
        ts=open(ts_file,'w')
        ts.write(timestamp)
        ts.close()


def download_taxdb(data_dir):

    """
    Download taxonomy database from NCBI (if newer than local version)
    and extracts names.dmp and nodes.dmp files

    arguments:
    none

    returns:
    none
    """

    download_file(data_dir, 'ftp.ncbi.nih.gov', '/pub/taxonomy', 'taxdump.tar.gz')

    if os.path.isfile('{}/{}'.format(data_dir,'taxdump.tar.gz')):
        tar=tarfile.open('{}/taxdump.tar.gz'.format(data_dir),'r:gz')
        tar.extract('nodes.dmp',path=data_dir)
        tar.extract('names.dmp',path=data_dir)
        os.remove('{}/taxdump.tar.gz'.format(data_dir))

    return()

def parse_taxonomy(data_dir):

    """
    Parse the nodes.dmp and names.dmp files to produce a dict containing
    the name, rank and parent node id of each node in the database

    arguments:
    data_dir -- path to directory containing taxonomy files

    returns:
    nodes -- dict of dicts
    """

    # The nodes.dmp file is delimited by '\t|\t' and contains the following
    # fields (from https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump_readme.txt)

    #tax_id                 -- the id of node associated with this name
    #name_txt               -- name itself
    #unique name                -- the unique variant of this name if name not unique
    #name class             -- (synonym, common name, ...)

    # The species names are found in the 'scientific name' lines

    names = dict()
    nodes = dict()

    for line in open(data_dir + '/names.dmp','r'):
        fields = line.rstrip("\n").split("\t|\t")
        if 'scientific name' in fields[3]:
            names[int(fields[0])] = fields[1]
    
    # The nodes.dmp file is delimited by '\t|\t' and contains the following
    # fields (from https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump_readme.txt)

    # tax_id -- node id in GenBank taxonomy database
    # parent tax_id             -- parent node id in GenBank taxonomy database
    # rank                  -- rank of this node (superkingdom, kingdom, ...)
    # embl code             -- locus-name prefix; not unique
    # division id               -- see division.dmp file
    # ....

    for line in open(data_dir + '/nodes.dmp', 'r'):
        fields = line.rstrip("\n").split("\t|\t")
        node_data = {
            'parent': int(fields[1]),
            'rank': fields[2],
            'name': names[int(fields[0])]
        }
        nodes[int(fields[0])] = node_data
    
    f=lambda x: x.replace('\t|','')
    tax_names=pd.read_csv('{}/names.dmp'.format(data_dir), sep="\t\|\t",header=0,engine='python',converters={3:f},
                      names=['taxid','name','dontknow','type'])
    tax_names=tax_names.loc[tax_names['type']=='scientific name']
    tax_names=tax_names[['taxid','name']]

    return(nodes,tax_names)

def get_tax_info_by_id(taxnodes,info,id):
    id=int(id)
    node=taxnodes[id]
    info[node['rank']]=node['name']
    if node['parent']==1:
        return(info)
   
    info=get_tax_info_by_id(taxnodes,info,node['parent'])
    return(info)

def get_taxnodes(data_dir):
    if not os.path.exists('{}/taxnodes.pkl'.format(data_dir)):
        print('Parsing taxnodes')
        name2taxid={}
        taxnodes,names=parse_taxonomy(data_dir)
        for key in taxnodes.keys():
            name2taxid[taxnodes[key]['name']]=key
        file= open('{}/taxnodes.pkl'.format(data_dir),'wb') 
        pickle.dump(taxnodes,file)
        file.close()
        file = open('{}/name2taxid.pkl'.format(data_dir),'wb') 
        pickle.dump(name2taxid, file)
        file.close()
    else:
        print('Loading existing taxnodes...')
        file = open('{}/taxnodes.pkl'.format(data_dir),'rb')
        taxnodes=pickle.load(file)
        file.close()
        file = open('{}/name2taxid.pkl'.format(data_dir),'rb')
        name2taxid=pickle.load(file)
        file.close()

    return(taxnodes,name2taxid)
    
def get_species_for_taxid(taxid, nodes):

    """
    Walk up taxonomy tree from defined taxid until the species node is
    found

    arguments:
    taxid -- taxonomy id
    nodes -- dict: parsed taxonomy nodes

    returns:
    node_data -- dict (success), string('missing/no_species') on error
    """

    if int(taxid) in nodes:
        node_data=nodes[int(taxid)]
        rank=node_data['rank']
        while node_data['rank'] != 'species':
            # seem to be some taxa which don't have a species defined...
            # catch these for reporting later...
            if node_data['parent'] == 1:
                return('no_species')
            node_data=nodes[node_data['parent']]
        return(node_data)

    else:
        return('missing')
