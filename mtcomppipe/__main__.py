import requests 
import re
import pandas as pd
import datetime
import random
import subprocess
import os
import sys
from colorama import Fore, Back, Style
from argparse import ArgumentParser, Namespace
from bs4 import BeautifulSoup

def print_error_message(message):
    """ Print error message in red
    parameters: Message to print
    return: None
    """
    print(Fore.RED + Style.BRIGHT + '[ERROR]', end=' ')
    print(message)
    print(Style.RESET_ALL)


def download_meta_data(id):
    """Download metadata by sample id

    parameters: Sample id
    return: Metadata dataframe
    """
    
    url ='https://amtdb.org/download_metadata'

    payload = {
        'csrfmiddlewaretoken': '0HqHvj9FNjUmQFUbV1Tnn37PAEqWTNaBGK9na0245R2GjRivATuSiiAxa2zYrjpI',
        'sample_ids': id
    }
    headers = {
        'User-Agent': 'Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/121.0.0.0 Safari/537.36',
        'Referer': 'https://amtdb.org/samples',
        'authority': 'amtdb.org',
        'method': 'POST',
        'path': '/download_metadata',
        'scheme': 'https',
        'Accept': 'text/html,application/xhtml+xml,application/xml;q=0.9,image/avif,image/webp,image/apng,*/*;q=0.8,application/signed-exchange;v=b3;q=0.7',
        'Accept-Encoding': 'gzip, deflate, br',
        'Accept-Language': 'en-US,en;q=0.9',
        'Cache-Control': 'max-age=0',
        'Content-Length': '106',
        'Content-Type': 'application/x-www-form-urlencoded',
        'Cookie': 'csrftoken=skvkjSVCb78sLFWnwa8Wp5rjkMp8PVKu8ne0YzO1tFgMeRkHb2JrkkU1UayanrZB; cookieconsent_status=dismiss; _gid=GA1.2.1056843561.1707838268; _ga=GA1.2.1007400171.1704845570; _ga_X26YEGG780=GS1.1.1707851898.14.0.1707851898.0.0.0',
        'Origin': 'https://amtdb.org',
        'Referer': 'https://amtdb.org/samples',
        'Sec-Ch-Ua': '"Not A(Brand";v="99", "Google Chrome";v="121", "Chromium";v="121"',
        'Sec-Ch-Ua-Mobile': '?0',
        'Sec-Ch-Ua-Platform': "Windows",
        'Sec-Fetch-Dest': 'document',
        'Sec-Fetch-Mode': 'navigate',
        'Sec-Fetch-Site': 'same-origin',
        'Sec-Fetch-User': '?1',
        'Upgrade-Insecure-Requests': '1',
    }

    r = requests.post(url, headers=headers, data=payload)
    raw_data = r.content.decode('utf-8')

    lines = [x.split(',') for x in raw_data.split('\n')]
    if len(lines) > 2:
        lines = lines[:2]

    df = pd.DataFrame(lines)
    df.columns = df.iloc[0]
    df = df[1:]
    df = df.reset_index(drop=True)

    return df    
 
def check_fastq_link(df):
    """ Check if project's metadata has a valid url to fastq
    parameters: Metadata dataframe
    return: True if it has a valid link, False otherwise
    """ 
    for value in df.iloc[0]:
        ebi_pattern = re.compile("^\"https?://www.ebi.ac.uk/*")
        ncbi_pattern = re.compile("^\"https?://www.ncbi.nlm.nih.gov/*")
        
        if ebi_pattern.match(value):
            url = value.replace('\"', '')
            alt_id = url.split('/')[-1]
            return True
        elif ncbi_pattern.match(value): 
            url = value.replace('\"', '')
            return False
    return False
  
def get_fastq_link(df):
    """ Look for url to fastq in metadata
    parameters: Metadata dataframe
    return: Url to fastq
    """  

    url = ''
    for value in df.iloc[0]:
        ebi_pattern = re.compile("^\"https?://www.ebi.ac.uk/*")
        if ebi_pattern.match(value):
            url = value.replace('\"', '')
    return url        

def get_fastq(alt_id, fasta_name, output_dir): 
    """ Download table of fastq associated with a project alternative id,
    and download those fastqs
    parameters: alternative id, fasta name and output dir
    return: total count of fastq, total size in bytes
    """

    # Request fastq page 
    url = f'https://www.ebi.ac.uk/ena/portal/api/filereport?accession={alt_id}&result=read_run&fields=run_accession,fastq_ftp,fastq_md5,fastq_bytes'
    page = requests.get(url) 

    # If can't access domain, abort and return None
    if page.status_code != 200:
        message = 'Status code ' + str(page.status_code) + ' while trying to connect to ' + url
        print_error_message(message)
        return None, None

    
    raw_data = page.content.decode('utf-8')
    lines = [x.split('\t') for x in raw_data.split('\n')]

    df = pd.DataFrame(lines)
    df.columns = df.iloc[0]
    df = df[1:]
    df = df.reset_index(drop=True)
    
    total_bytes = 0
    ras_count = 0
    for i in range(len(df)):
        row = df.iloc[i]
        run_accession = df.iloc[i]['run_accession']
        fastq_bytes = df.iloc[i]['fastq_bytes']
        
        if run_accession != None and fastq_bytes != None:
            if ';' in fastq_bytes:
                fastq_bytes = fastq_bytes.split(';')
                fastq_bytes = str(sum([int(i) in fastq_bytes]))
            try:
                # Export sra-tools to path
                env = os.environ.copy()
                env['PATH'] = env['PATH'] + ':' + os.path.abspath(os.path.join(os.path.dirname(__file__), 'sratoolkit.3.0.10-ubuntu64/bin/'))
                
                # Download fastq using parallel-fastq-dump
                subprocess.run(['parallel-fastq-dump', '--sra-id', run_accession, '--threads', '4', '--outdir', output_dir,'--split-files', '--gzip'], env=env)
                total_bytes += int(fastq_bytes)
                ras_count += 1
            except SystemExit:
                print_error_message('Could not download file with run acession:', run_accession)
                continue

            # fastq_name = run_accession
            # fastq_name += '.fastq.gz'
            
            # subprocess.run(['classpipe', '--refDNA', os.path.join(output_dir, fasta_name), '--aDNA1', fastq_name,  '--output', output_dir_fastq, '--saveBam'])
            # subprocess.run(['rm', '-rf', fastq_name])
    
    return ras_count, total_bytes

def get_info_from_amtdb():
    """ Download table of samples from amtdb, 
    and filter rows that contain fastq links
    parameters: None
    return: sample id list, link to fasta list, metadata map
    """

    url = 'https://amtdb.org/samples'
    page = requests.get(url) 
    
    if page.status_code != 200:
        message = 'Status code ' + str(page.status_code) + ' while trying to connect to ' + url
        print_error_message(message)
        exit(1)

    soup = BeautifulSoup(page.content, 'html.parser')
    table = soup.find('table', {'id': 'table-samples'})
    
    # Get table rows and remove header
    rows = table.findAll('tr')
    rows = rows[1:]
    
    # Sample Id
    ids = []
    # Sample link to fasta
    links = []

    for row in rows:
        # Get '<a> tag' with fasta download link 
        aTag = row.find('a', {'class': 'a-control a-download_fasta text-decoration-none'})
        if aTag != None:
            # Mount link to download sample
            link = 'https://amtdb.org' + aTag['href']
            # Add id and link to lists
            links.append(link)
            ids.append(row['sample_id'])
    
    # Start random generator with timestamp
    indexes = []
    
    md_map = {}
    # Get list of random indexes
    for i in range(len(ids)):
        metadataDF = download_meta_data(ids[i])
        validFastqLink = check_fastq_link(metadataDF)
        
        if validFastqLink and ids[i] not in indexes:
            md_map[ids[i]] = metadataDF
            indexes.append(i)

    # Filter list by generated indexes
    ids = [ids[i] for i in indexes]
    links = [links[i] for i in indexes]
    
    return ids, links, md_map

def main():
    """ 
    """
    
    parser = ArgumentParser()

    # Set number of samples from command line
    parser.add_argument('-n', help='Number of references to download')
    parser.add_argument('--output', help='Set output destination [Default: .]', default='.')

    # Get arguments from command line
    args: Namespace = parser.parse_args()

    # Number of samples requested
    nSamples = int(args.n)

    # Get samples table from AMTDB website
    ids, links, md_map = get_info_from_amtdb()

    if nSamples > len(ids):
        nSamples = len(ids)
        
    # Create output dir if it doesn't exist
    output_dir = args.output
    if args.output == '.':
        output_dir = os.path.join(args.output, 'samples')
    subprocess.run(['rm', '-rf', output_dir])
    subprocess.run(['mkdir', output_dir])

    indexes = [i for i in range(len(ids))]
    random.shuffle(indexes)
    
    total_fastq_count = 0
    total_fastq_size = 0
    projects_count = 0
    i=0
    # Download samples to output directory
    while projects_count < nSamples:
        # Create output dir for files
        output_dir_sample = os.path.join(output_dir, ids[indexes[i]])
        subprocess.run(['rm', '-rf', output_dir_sample])
        subprocess.run(['mkdir', output_dir_sample])
        
        # Download FASTA
        subprocess.run(['wget', '-P', output_dir_sample, links[indexes[i]]])
        
        # Save metadata to file
        metadataDF = md_map[ids[indexes[i]]]
        # print(md_map)
        # exit(0)
        # metadataDF = download_meta_data(ids[indexes[i]])
        metadataFile = open(os.path.join(output_dir_sample, './metadata.txt'), 'w')
        metadataFile.write(metadataDF.to_string())
        metadataFile.close()
        
        # Download FASTQs
        data_link = get_fastq_link(metadataDF)
        
        alt_id = data_link.split('/')[-1]
        alt_id = alt_id.split('?')[0]
        fasta_name = links[indexes[i]].split('/')[-1]
        fastq_count, size_bytes = get_fastq(alt_id, fasta_name, output_dir_sample)
        
        if fastq_count and size_bytes:
            print('=====================')
            print(alt_id)
            print('\tFASTQ count:\t', fastq_count)
            print('\tFASTQ size (Gb):\t', size_bytes/(1024**3))
            total_fastq_count += fastq_count
            total_fastq_size += size_bytes
            projects_count += 1
            i += 1
    
    print('=====================')
    print('\tTotal FASTQ count:\t', total_fastq_count)
    print('\tTotal FASTQ size (Gb):\t', total_fastq_size/(1024**3))

    # Mount pandas dataframe and filter by samples in ids list 
    # df = pd.read_html(str(table))[0]
    df = df[df['Name'].isin(ids)]

    # Get current year 
    # year = datetime.date.today().year
    # Function to calculate estimated average age based on range (Columns 'Year from' and 'Year to')
    # fun = lambda x: year + (abs(x['Year from']) + abs(x['Year to']))/2 if x['Year from'] < 0 and x['Year to'] < 0 else year - (x['Year from'] + x['Year to'])/2

    # Apply lambda function to dataframe and create new column from results
    # df['Age'] = df.apply(fun, axis=1)

    # print(df.columns)
    # print(df[['Name', 'Year from', 'Year to', 'Age']])

    # Save metadata to file
    
if __name__=='__main__':
    main()