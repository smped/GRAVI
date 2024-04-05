## x is config['genome']
def get_gtf_url(x):
    bld = x['build'].lower()
    sp_map = {
        'grch37': "human", 'grch38': "human",
        'hg19': "human", 'hg38': "human",
        'grcm39': 'mouse', 'grcm38': 'mouse',
        'mm10': 'mouse', 'mm39': 'mouse'
    }
    sp = sp_map[bld]
    gc_vers = re.findall('[M0-9]+$', x['version'])[0]
    subdir = ''
    extra = ''
    if (bld == 'grch37') | (bld == 'hg19'):
        subdir = '/GRCh37_mapping'
        extra = 'lift37'
    url = list(urllib.parse.urlparse("https://ftp.ebi.ac.uk"))
    url[2] = "pub/databases/gencode/Gencode_" + sp + "/release_" + gc_vers + subdir + '/' + gtf
    return(urllib.parse.urlunparse(url))


rule download_gtf:
    output: gtf
    params:
        url = get_gtf_url(config['genome'])
    threads: 1
    localrule: True
    log: os.path.join(log_path, "downloads", "gtf.log")
    shell:
        """
        curl --fail-early \
            --output {output} \
            {params.url} 2> {log}
        """

rule download_blacklist:
    output: blacklist
    params:
        url = "https://github.com/Boyle-Lab/Blacklist/blob/master/lists/" + os.path.basename(blacklist)
    threads: 1
    localrule: True
    log: os.path.join(log_path, "downloads", "blacklist.log")
    shell:
        """
        curl --fail-early \
            --output {output} \
            {params.url} 2> {log}
        """