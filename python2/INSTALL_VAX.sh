#!/bin/sh

#download VAX data and import into mysql database

#VAX requires a MySQL client, access to a MySQL server, and a full local installation of the Ensembl API.

#This installer requires ncftpget and wget

#Some plugins, such as Alignment and Protein, require a local Ensembl MySQL database for best performance.

#Create PERL5LIB references to the ensembl and bioperl-live 1.2.3 modules,
#a reference to Bio::Matrix::IO from a current BioPerl download,
#and a reference to SWISS::Knife from http://sourceforge.net/projects/swissknife/files/latest/download.

#This installer uses passwords in the command line for database creation and access.
#The MySQL 5.6 client generates a warning, which may be ignored depending on your security concerns.
#For more information and possible workarounds, see http://bugs.mysql.com/bug.php?id=66546
#and http://akrabat.com/software/password-less-command-line-scripts-with-mysql-5-6/

#get path and parent directory of this script
#http://hintsforums.macworld.com/archive/index.php/t-73839.html
# if the call came via symlink, then use its target instead:
arg=$0; [[ -L $0 ]] && arg=$(stat -L -c'%n' "$0");
script_path=$(2>/dev/null cd "${arg%/*}" >&2; echo "`pwd -P`/${arg##*/}");
script_dir=$(dirname "$script_path");

HERE=${script_dir};

wd=$(pwd);

echo "VAX installer: Plugins and database for use with the Ensembl Variant Effect Predictor.";
echo "Runtime about 1 hour with a fast internet connection (>1 day if -S option used)";
echo;
echo "For help: $0 -?";
echo;

function pause(){
    t=$1; shift;
   read -t $t -p "$*"
}

function log(){
    echo $(date "+%Y%m%d %H:%M:%S: ") "$*";
}

#requirements
for r in {mysql,ncftpget,wget,curl}; do
    if [ ! -x "$(which $r)" ]; then
        echo "$r executable required. Aborting." exit 1;
    fi
done

#defaults
HOST=localhost;
PORT=3306;
USER=;
PASSWORD=;
DATABASE=vax;
VAX_USER=vax;
VAX_PW=vax;
VEP_DIR=${HOME}/.vep;
CONFIG_FILE=${VEP_DIR}/vep.ini;
PLUGINS_DIR=${VEP_DIR}/Plugins;
ENSEMBL_VERSION=0; #75
VEP_DESTDIR=${VEP_DIR}/VEP;
VEP_CACHEDIR=${VEP_DIR};
INSTALL_VEP=0;
INSTALL_ENSEMBL_API=0;
INSTALL_ENSEMBL_DATABASES=0
BIOPERL_VERSION="release-1-6-923";
VEP_AUTO=;
VEP_SPECIES="--SPECIES homo_sapiens";
DOWNLOAD_DIR=${HOME}/vax_downloads;
DOWNLOAD=1;
INSTALL_DB=1;
INSTALL_VAX_PLUGINS=1;
INSTALL_VEP_PLUGINS=0;
DOWNLOAD_FATHMM=0;
DOWNLOAD_DBNSFP_VERSION=0; #2.1
DBNSFP_DIR=${VEP_DIR}/dbnsfp;
HGMD_USER=hgmd;
HGMD_PW=hgmd;

help_file="Usage: $0 [options]
Options:
    [-h] <msyql_host> (default: ${HOST})
    [-P] <mysql_port> (default: ${PORT})
    [-u] <mysql_install_user> (default: enter at prompt)
    [-p] <mysql_install_password> (default: enter silently at prompt)
    [-d] <mysql_database> (default: ${DATABASE}) EXISTING DATABASE WILL BE OVERWRITTEN!
    [-U] <mysql user for runtime access to database> (default: ${VAX_USER})
    [-W] <mysql password for runtime access> (default: ${VAX_PW})
    [-v] </path/to/.vep directory> (default: ${VEP_DIR})
    [-D] </path/to/vax_downloads directory> (default: ${DOWNLOAD_DIR})
    [-X] (install VEP_plugins downloaded from git; default: ${INSTALL_VEP_PLUGINS} = do not download/install; beta test-not a supported VAX function)
    [-e] <version of Ensembl API/VEP/cache to install> (default: ${ENSEMBL_VERSION} = do not download/install)
    [-E] <VEP destination directory for Ensembl API and Variant Effect Predictor> (default: ${VEP_DESTDIR})
    [-I] (run VEP installer to download and install subset of Ensembl API, cache, and fasta; requires -e; default: ${INSTALL_VEP} = do not download/install; beta test-not a supported VAX function)
    [-i] (download the entire Ensembl Perl API, and a current version of BioPerl; requires -e; default: ${INSTALL_ENSEMBL_API} = do not download; beta test-not a supported VAX function. WARNING: The full Ensembl API and a current version of BioPerl may be required by some VAX plugins.)
    [-J] (download and install ensembl_ and homo_asapiens_ databases; default: ${INSTALL_ENSEMBL_DATABASES} = do ot install.)
    [-b] <version of BioPerl to download> (default: ${BIOPERL_VERSION}; requires -e and -i)
    [-a] (VEP --AUTO parameters; Run VEP installer without user prompts. Use a (API), c (cache), f (FASTA) to specify parts to install e.g., -a ac for API and cache; -a for cache only is recommended for VAX unless cashe is already installed; default: None)
    [-s] (VEP --SPECIES parameters; Comma-separated list of species to install when using --AUTO; default: homo_sapiens)
    [-F] (download and install FATHMM database; default: ${DOWNLOAD_FATHMM} = do not download/install; beta test-not a supported VAX function)
    [-S] <version of dbNSFP to download> (large ~3Gb; default: ${DOWNLOAD_DBNSFP_VERSION} = do not download/install;  beta test-not a supported VAX function)
    [-g] <mysql user for runtime access to hgmd_pro> (default: ${HGMD_USER})
    [-G] <mysql password for runtime access to hgmd_pro> (default: ${HGMD_PW})
    [-x] (do not download any files; assumes all files exist in ${DOWNLOAD_DIR} and have been pre-processed as needed; default: ${DOWNLOAD} = do download/install)
    [-y] (do not install MySQL database ${DATABASE}; default: ${INSTALL_DB} = do download/install)
    [-z] (do not create ${VEP_DIR}, ${CONFIG_FILE}, and ${PLUGINS_DIR}; default: ${INSTALL_VAX_PLUGINS} = do download/install)
    [-?] (print this usage info and exit)"

usage() {
    printf "%b\n" "${help_file}";
    1>&2; exit 1;
}

while getopts ":h:P:u:p:d:U:W:v:D:Xe:E:IiJba:s:FS:g:G:xyz" o; do
    case "${o}" in
        h)  HOST="${OPTARG}" ;;
        P)  PORT="${OPTARG}" ;;
        u)  USER="${OPTARG}" ;;
        p)  PASSWORD="${OPTARG}" ;;
        d)  DATABASE="${OPTARG}" ;;
        U)  VAX_USER="${OPTARG}" ;;
        W)  VAX_PW="${OPTARG}" ;;
        v)  VEP_DIR="${OPTARG}" ;;
        D)  DOWNLOAD_DIR="${OPTARG}" ;;
        X)  INSTALL_VEP_PLUGINS=1 ;;
        e)  ENSEMBL_VERSION="${OPTARG}" ;;
        E)  VEP_DESTDIR="${OPTARG}" ;;
        I)  INSTALL_VEP=1 ;;
        i)  INSTALL_ENSEMBL_API=1 ;;
        J)  INSTALL_ENSEMBL_DATABASES=1 ;;
        b)  BIOPERL_VERSION="${OPTARG}" ;;
        a)  VEP_AUTO="--AUTO ${OPTARG}" ;;
        s)  VEP_SPECIES="--SPECIES ${OPTARG}" ;;
        F)  DOWNLOAD_FATHMM=1 ;;
        S)  DOWNLOAD_DBNSFP_VERSION="${OPTARG}" ;;
        g)  HGMD_USER="${OPTARG}" ;;
        G)  HGMD_PW="${OPTARG}" ;;
        x)  DOWNLOAD=0 ;;
        y)  INSTALL_DB=0 ;;
        z)  INSTALL_VAX_PLUGINS=0 ;;
        *)  usage; exit ;;
    esac
done
shift $((OPTIND-1))

if [[ ${INSTALL_DB} == 1 ]]; then
    if [ -z "${USER}" ]; then
        read -p "Enter MySQL user name for database installation (or run using the -u option): " -r
        #echo    # (optional) move to a new line
        if [ -z ${REPLY} ]
        then
            echo;
            log "Installer aborted."
            exit 1;
        else
            USER=${REPLY};
        fi
    fi
    
    if [ -z "${PASSWORD}" ]; then
        read -p "Enter MySQL password for database installation (or run using the -p option): " -s -r
        #echo    # (optional) move to a new line
        if [ -z ${REPLY} ]
        then
            echo;
            log "Installer aborted."
            exit 1;
        else
            PASSWORD=${REPLY};
        fi
    fi

    if [ -z "${USER}" ] || [ -z "${PASSWORD}" ]; then
        echo "MySQL user (-u) and password (-p) parameters are required."
        usage;
        exit;
    fi
fi

if [[ "${VEP_DIR}" != "${HOME}/.vep" ]]; then
    CONFIG_FILE=${VEP_DIR}/vep.ini;
    PLUGINS_DIR=${VEP_DIR}/Plugins;
    VEP_CACHEDIR=${VEP_DIR};
    DBNSFP_DIR=${VEP_DIR}/dbnsfp;
fi

#MySQL commands
MYSQL="$(which mysql) -h ${HOST} -P ${PORT} --user=${USER} -p${PASSWORD}";
MYSQLIMPORT="$(which mysqlimport) -h ${HOST} -P ${PORT} --user=${USER} -p${PASSWORD}";

#directories
ENSEMBL_API_DIR=${VEP_DESTDIR}/ensembl-${ENSEMBL_VERSION};
BIOPERL_DIR=${VEP_DESTDIR}/bioperl-${BIOPERL_VERSION};

echo; echo;
#print config
cat << EOF
Parameters:
    HOST (-h) = ${HOST}
    PORT (-P UPPERCASE) = ${PORT}
    USER (-u) = ${USER}
    PASSWORD (-p lowercase) = ********
    DATABASE (-d) = ${DATABASE}
    VAX_USER (-U) = ${VAX_USER}
    VAX_PW (-W) = ${VAX_PW}
    VEP_DIR (-v) = ${VEP_DIR}
    CONFIG_FILE = ${CONFIG_FILE}
    PLUGINS_DIR = ${PLUGINS_DIR}
    ENSEMBL_VERSION (-e) = ${ENSEMBL_VERSION}
    VEP_DESTDIR (-E) = ${VEP_DESTDIR}
    VEP_CACHEDIR = ${VEP_CACHEDIR}
    INSTALL_VEP (-I) = ${INSTALL_VEP}
    INSTALL_ENSEMBL_API (-i) = ${INSTALL_ENSEMBL_API}
    BIOPERL_VERSION (-b) = ${BIOPERL_VERSION}
    VEP_AUTO (-a) = ${VEP_AUTO}
    VEP_SPECIES (-s) = ${VEP_SPECIES}
    DOWNLOAD_DIR (-D) = ${DOWNLOAD_DIR}
    DOWNLOAD (-x) = ${DOWNLOAD}
    INSTALL_DB (-y) = ${INSTALL_DB};
    INSTALL_VAX_PLUGINS (-z) = ${INSTALL_VAX_PLUGINS}
    INSTALL_VEP_PLUGINS (-X) = ${INSTALL_VEP_PLUGINS}
    DOWNLOAD_FATHMM (-F) = ${DOWNLOAD_FATHMM}
    DOWNLOAD_DBNSFP_VERSION (-S) = ${DOWNLOAD_DBNSFP_VERSION}
    DBNSFP_DIR = ${DBNSFP_DIR}
    HGMD_USER (-g) = ${HGMD_USER}
    HGMD_PW (-G) = ${HGMD_PW}
EOF

echo;
if [[ ${DOWNLOAD} != 1 ]]; then
    echo "  Remote data will not be downloaded.";
    echo "    All files must exist in ${DOWNLOAD_DIR} from a prior run.";
else
    echo "  Remote data will be downloaded.";    
fi

echo;
if [[ ${INSTALL_DB} != 1 ]]; then
    echo "  MySql database for VAX plugins will not be installed";
else
    echo "  MySql database ${DATABASE} for VAX plugins will be installed.";
    echo "    ANY EXISTING ${DATABASE} DATABASE WILL BE DELETED AND REPLACED!";
    echo "    If you want to download data but not install the database, use the -y option."
fi

echo;
if [[ ${INSTALL_VAX_PLUGINS} != 1 ]]; then
    echo "  ${VEP_DIR}, ${CONFIG_FILE}, and ${PLUGINS_DIR} will not be installed.";
else
    cat << EOF
  The following objects will be created if they do not exist:
    ${VEP_DIR}
    ${CONFIG_FILE}
    ${PLUGINS_DIR}
  These objects will be backed up:
    ${CONFIG_FILE}
    ${PLUGINS_DIR}
  If you want to download data and/or install the database
  but not install these objects, use the -z option.
EOF
fi

echo;
if [[ ${INSTALL_VEP_PLUGINS} == 1 ]]; then
    echo "  VEP_plugins downloaded from git will be installed (this is a beta test function not VAX-supported).";
    if [[ ${DOWNLOAD_FATHMM} == 1 ]]; then
        echo "  FATHMM plugin database will be downloaded and installed on the mysql server (this is a beta test function not VAX-supported).";
    fi
    if [[ ${DOWNLOAD_DBNSFP_VERSION} != 0 ]]; then
        echo;
        echo "  dbNSFP plugin data will be downloaded and installed (this is a beta test function not VAX-supported).";
    fi
fi

echo;
if [[ ${ENSEMBL_VERSION} != 0 ]]; then
    echo;
    echo "The Variant Effect Predictor script will be downloaded as:
    ${VEP_DESTDIR}/variant_effect_predictor/variant_effect_predictor.pl"
    if [[ ${INSTALL_VEP} == 1 ]]; then
        echo;
        echo "The VEP installer will be executed with the following parameters:";
        echo "  --DESTDIR ${VEP_DESTDIR} --CACHEDIR ${VEP_CACHEDIR} ${VEP_SPECIES} ${VEP_AUTO}"
    else
        echo "The VEP installer will not be executed.";
    fi
    if [[ ${INSTALL_ENSEMBL_API} == 1 ]]; then
        echo;
        echo "Version ${ENSEMBL_VERSION} of the Ensembl Perl API will be downloaded to ${ENSEMBL_API_DIR}.";
        echo "BioPerl ${BIOPERL_VERSION} will be downloaded to ${BIOPERL_DIR}.";        
    else
        echo "The Ensembl Perl API and BioPerl will not be downloaded.";
    fi
fi
echo;
if [[ ${INSTALL_DB} == 1 ]]; then
    #test MySQL server connection
    log "testing MySQL connection.";
    echo;
    cat << EOF
    MySQL client 5.6 will issue password warnings.
    In most situations these may safely be ignored.
    See the README for more information.
EOF

    echo;
    if ${MYSQL} -e ";" ; then
        log "MySQL connection OK.";
    else
        log "MySQL connection failed.";
        echo "  Verify parameters for HOST(-h), PORT(-P), USER(-u), and PASSWORD(-p).";
        log "Installer aborted."
        exit 1;
    fi
fi

#user must approve parameters
echo;
read -p "Continue? y [n] " -n 1 -r
#echo    # (optional) move to a new line
if [[ ! $REPLY =~ ^[Yy]$ ]]
then
    echo;
    log "Installer aborted."
    exit 1;
fi

################################################################################
#                                    debug                                     #
################################################################################
#exit;

echo;
log "Starting VAX installation."

################################################################################
#                             create directories                               #
################################################################################

echo;
#directories
log "creating directories:";
echo "  ${DOWNLOAD_DIR}";
mkdir -p ${DOWNLOAD_DIR};
if [[ ${INSTALL_VAX_PLUGINS} == 1 ]]; then
    echo "  ${VEP_DIR}";
    mkdir -p ${VEP_DIR};
    echo "  ${PLUGINS_DIR}";
    mkdir -p ${PLUGINS_DIR};
fi
if [[ ${ENSEMBL_VERSION} != 0 ]]; then
    echo "  ${VEP_DESTDIR}";
    mkdir -p ${VEP_DESTDIR};
    echo "  ${VEP_CACHEDIR}";
    mkdir -p ${VEP_CACHEDIR};
fi
if [[ ${DOWNLOAD_DBNSFP_VERSION} != 0 ]]; then
    echo "  ${DBNSFP_DIR}";
    mkdir -p ${DBNSFP_DIR};
fi

################################################################################
#                      drop and re-create VAX database                         #
################################################################################

if [[ ${INSTALL_DB} == 1 ]]; then
    echo;
    log "dropping database ${DATABASE}";
    #fail-safe
    read -p "Continue? y [n] " -n 1 -r
    if [[ ! $REPLY =~ ^[Yy]$ ]]
    then
        echo;
        log "Installer aborted."
        exit 1;
    fi
    echo;
    ${MYSQL} -e "DROP DATABASE IF EXISTS ${DATABASE};";
    log "creating database ${DATABASE}";
    ${MYSQL} -e  "CREATE DATABASE ${DATABASE};";
    log "granting all permissions on ${DATABASE} to ${VAX_USER} from anywhere with password ${VAX_PW}";
    ${MYSQL} -e  "GRANT ALL ON ${DATABASE}."'*'" TO '${VAX_USER}'@'%' IDENTIFIED BY '${VAX_PW}';";
fi

################################################################################
#                               Ensembl API & VEP                              #
################################################################################

BAK_SUFFIX="_"$(date +%Y%m%d%H%M%S);
if [[ ${ENSEMBL_VERSION} != 0 ]]; then
    echo;
    log "Downloading and installing Variant Effect Predictor to ${VEP_DESTDIR}";
    log "The Variant Effect Predictor executable will be installed at:";
    log "  ${VEP_DESTDIR}/variant_effect_predictor/variant_effect_predictor.pl"
    if [ -d ${VEP_DESTDIR} ]; then
        log "Backing up existing ${VEP_DESTDIR} to ${VEP_DESTDIR}${BAK_SUFFIX}";
        cp -a ${VEP_DESTDIR} ${VEP_DESTDIR}${BAK_SUFFIX};
    fi
    mkdir -p ${VEP_DESTDIR};
    mkdir -p ${VEP_CACHEDIR};
    this_dir=$(pwd); cd ${VEP_DESTDIR};
    curl "http://cvs.sanger.ac.uk/cgi-bin/viewvc.cgi/ensembl-tools/scripts/variant_effect_predictor.tar.gz?view=tar&root=ensembl&pathrev=branch-ensembl-${ENSEMBL_VERSION} -o -" | tar xvz -C ${VEP_DESTDIR};
    if [[ -d ${VEP_DESTDIR}/variant_effect_predictor && -e ${VEP_DESTDIR}/variant_effect_predictor/INSTALL.pl ]]; then
        if [[ ${INSTALL_VEP} == 1 ]]; then
            echo "WARNING: during VEP installation, consider downloading cache but not the Ensembl/BioPerl API subset and the fasta file(s), which may not work well with some VAX plugins. You can use the -i option to get full versions of the Perl modules."
            cd ${VEP_DESTDIR}/variant_effect_predictor;
            perl ${VEP_DESTDIR}/variant_effect_predictor/INSTALL.pl --DESTDIR ${VEP_DESTDIR} --CACHEDIR ${VEP_CACHEDIR} ${VEP_SPECIES} ${VEP_AUTO};
        fi
    else
        log "failed to download VEP installer ${VEP_DESTDIR}/variant_effect_predictor/INSTALL.pl."
        echo "maybe you should try installing VEP manually";
        exit 0;
    fi
    cd ${this_dir};
    
    #full Ensembl Perl API and BioPerl
    if [[ ${INSTALL_ENSEMBL_API} == 1 ]]; then
        echo;
        log "downloading Ensembl Perl API version ${ENSEMBL_VERSION} to ${ENSEMBL_API_DIR}";
        this_dir=$(pwd);
        if [ -d ${ENSEMBL_API_DIR} ]; then
            log "deleting existing ${ENSEMBL_API_DIR}."
            rm -rf ${ENSEMBL_API_DIR};
        fi
        mkdir -p ${ENSEMBL_API_DIR};
        cd ${VEP_DESTDIR};
        rm ${VEP_DESTDIR}/ensembl;
        ln -s ensembl-${ENSEMBL_VERSION} ensembl;
        cd ${ENSEMBL_API_DIR};
        cvs -d :pserver:cvsuser@cvs.sanger.ac.uk:/cvsroot/ensembl checkout -r branch-ensembl-${ENSEMBL_VERSION} ensembl-api;
        cd ${this_dir};
        
        echo;
        log "downloading BioPerl ${BIOPERL_VERSION} to ${BIOPERL_DIR}";
        this_dir=$(pwd);
        if [ -d ${BIOPERL_DIR} ]; then
            log "deleting existing ${BIOPERL_DIR}."
            rm -rf ${BIOPERL_DIR};
        fi
        mkdir -p ${BIOPERL_DIR};
        cd ${VEP_DESTDIR};
        rm ${VEP_DESTDIR}/bioperl-live;
        ln -s bioperl-${BIOPERL_VERSION} bioperl-live;
        cd ${BIOPERL_DIR};
        git clone git://github.com/bioperl/bioperl-live.git ${BIOPERL_DIR};
        cd ${BIOPERL_DIR};
        #get list of BioPerl versions
        #git tag -l
        git checkout ${BIOPERL_VERSION};
        cd ${this_dir};
    fi
fi

if [[ ${INSTALL_VAX_PLUGINS} == 1 ]]; then
    if [ -d ${PLUGINS_DIR} ]; then
        log "Backing up existing ${PLUGINS_DIR} to ${PLUGINS_DIR}${BAK_SUFFIX}";
        cp -a ${PLUGINS_DIR} ${PLUGINS_DIR}${BAK_SUFFIX};
    fi
    if [ -e ${CONFIG_FILE} ]; then
        log "Backing up existing ${CONFIG_FILE} to ${CONFIG_FILE}${BAK_SUFFIX}";
        mv ${CONFIG_FILE} ${CONFIG_FILE}${BAK_SUFFIX};
    fi

    #copy Plugins to user's Plugin directory
    log "installing vax plugins and config file:"
    echo "  ${PLUGINS_DIR}";
    cp -a ${HERE}/Plugins/* ${PLUGINS_DIR}/;
    
    #write vax_vep.ini to user's .vep
    echo "  ${CONFIG_FILE}";
    cat > ${CONFIG_FILE} << EOF
#VAX plugins
#Edit this file for your needs, add the plugin options to your own vep.ini,
#or include the plugins (--plugin PluginName) on the VEP command line
#${PLUGINS_DIR}/VAX.pm must be referenced in PERL5LIB
#for example by adding the following line (without the #) to ${HOME}/.bashrc 
#export PERL5LIB=${HOME}/.vep/Plugins:${PERL5LIB}
#
#The vw plugin is required by other plugins that use the MySQL database:
#  GeneIDs,HPA,KEGG,MousePhenotypes,OMIM,UniProt
#
#change host, port, user, password to your local Ensembl database
#or omit for remote Ensembl (not recommended) or cache only
database=1
host=${HOST}
port=${PORT}
user=ensembl
password=ensembl
cache=1
#
#plugin vw must preceed any plugins that use it
#
plugin=vw,${HOST},${PORT},${VAX_USER},${VAX_PW},mysql,${DATABASE}
plugin=Alignment
plugin=Consequences
plugin=DiseasesPhenotypes
plugin=GeneIDs
##HGMD requires commercial hgmd_pro database and stored procedures coord2hgmd and gene2hgmd_disease
#plugin=HGMD,${HOST},${PORT},${HGMD_USER},${HGMD_PW},mysql,hgmd_pro
plugin=HPA
plugin=KEGG
plugin=MousePhenotypes
plugin=OMIM
plugin=Phenotypes
plugin=Protein
plugin=UniProt
plugin=VCFCols
#typical options
verbose=1
force_overwrite=1
no_progress=1
##compress='gunzip -c' required for Mac OS X unless zcat installed
compress='gunzip -c'
#options for more output data
polyphen=b
sift=b
allele_number=1
biotype=1
canonical=1
ccds=1
check_alleles=1
check_existing=1
check_ref=1
#check_svs requires database
check_svs=1
core_type=core
domains=1
gene=1
gmaf=1
hgvs=1
#maf_1kg requires cache
maf_1kg=1
#maf_esp requires cache
maf_esp=1
#pubmed requires cache
pubmed=1
numbers=1
protein=1
regulatory=1
symbol=1
terms=so
xref_refseq=0
#optional downloaded VEP_plugins
#plugin=Carol
#plugin=Condel,${PLUGINS_DIR}/config/Condel/config,b,2
#plugin=Conservation
#plugin=GO
#plugin=Grantham
#FATHMM requires Python MySQLdb module
#plugin=FATHMM,'python ${PLUGINS_DIR}/fathmm/fathmm.py'
#plugin=dbNSFP,${DBNSFP_DIR}/dbNSFP.gz,Interpro_domain,SLR_test_statistic,SIFT_score_converted,LRT_score,LRT_score_converted,LRT_pred,MutationTaster_score,MutationTaster_score_converted,MutationTaster_pred,MutationAssessor_score,MutationAssessor_score_converted,MutationAssessor_pred,FATHMM_score_converted,FATHMM_pred,GERP++_NR,GERP++_RS,phyloP,29way_pi,29way_logOdds,LRT_Omega,1000Gp1_AC,1000Gp1_AF,1000Gp1_AFR_AC,1000Gp1_AFR_AF,1000Gp1_EUR_AC,1000Gp1_EUR_AF,1000Gp1_AMR_AC,1000Gp1_AMR_AF,1000Gp1_ASN_AC,1000Gp1_ASN_AF,ESP6500_AA_AF,ESP6500_EA_AF
#
#ExtraCols must follow any plugins that generate EXTRAS keys
plugin=ExtraCols
#
EOF

    echo;
    cat << EOF
You must include a reference to ${PLUGINS_DIR} in the PERL5LIB environment variable
as well as references to the Ensembl Perl API and  to BioPerl. For example, add the following to ${HOME}/.bashrc:
    export PERL5LIB="${BIOPERL_DIR}:\${PERL5LIB}"
    export PERL5LIB="${ENSEMBL_API_DIR}/ensembl-external/modules:\${PERL5LIB}"
    export PERL5LIB="${ENSEMBL_API_DIR}/ensembl-compara/modules:\${PERL5LIB}"
    export PERL5LIB="${ENSEMBL_API_DIR}/ensembl-functgenomics/modules:\${PERL5LIB}"
    export PERL5LIB="${ENSEMBL_API_DIR}/ensembl-variation/modules:\${PERL5LIB}"
    export PERL5LIB="${ENSEMBL_API_DIR}/ensembl/modules:\${PERL5LIB}"
EOF
    echo;

    if [[ ${CONFIG_FILE} == ${HOME}/.vep/vep.ini ]]; then
            echo "Variant Effect Predictor will use the installed config file by default."
    else
        echo "Add '--config ${CONFIG_FILE}' to the Variant Effect Predictor command line."
    fi
    echo;
    echo "If you change VAX host, port, user, password, or db name, make conforming changes to plugin=vw in the config file."
    echo;
fi

################################################################################
#             VEP_plugins examples from git (beta - not part of VAX)           #
################################################################################

if [[ ${INSTALL_VEP_PLUGINS} == 1 ]]; then
    #add code here to download and install additional plugins
    
    #download Ensembl VEP_plugins from git
    log "downloading and installing VEP_plugins"
    #delete any existing git directory for ensembl plugins
    ensembl_plugin_git_dir=${DOWNLOAD_DIR}/ensembl_plugin_git;
    if [[ -d ${ensembl_plugin_git_dir} ]]; then rm -rf ${ensembl_plugin_git_dir}; fi
    #create new git directory
    mkdir -p ${ensembl_plugin_git_dir};
    this_dir=$(pwd);
    cd ${ensembl_plugin_git_dir};
    #clone the remote get plugins directory
    git clone https://github.com/ensembl-variation/VEP_plugins.git;
    #copy the downloaded plugins to the VEP_plugins sirectory
    cp -a ${ensembl_plugin_git_dir}/VEP_plugins/* ${PLUGINS_DIR}/;
    
    #Condel plugin
    #put local path in Condel configuration file
    log "updating ${PLUGINS_DIR}/config/Condel/config/condel_SP.conf"
    #sed "s:path/to/config/Condel/:${PLUGINS_DIR}/config/Condel/:" ${PLUGINS_DIR}/config/Condel/config/condel_SP.conf > ${PLUGINS_DIR}/config/Condel/config/condel_SP.conf;
    cat > ${PLUGINS_DIR}/config/Condel/config/condel_SP.conf << EOF
condel.dir='${PLUGINS_DIR}/config/Condel/'

#------------------------------------------------------------------------------
cutoff.HumVar.sift='0.15'
cutoff.HumVar.polyphen='0.28'
cutoff.HUmVar.condel='0.46'
#------------------------------------------------------------------------------
max.HumVar.sift='1'
max.HumVar.polyphen='1'
EOF

    #dbNSFP plugin
    if [[ ${DOWNLOAD_DBNSFP_VERSION} != 0 ]]; then
        if [[ -e ${DBNSFP_DIR}/dbNSFP.gz && -e ${DBNSFP_DIR}/dbNSFP.gz.tbi ]]; then
            echo;
            echo "It looks like the dbSNFP data has been installed in ${DBNSFP_DIR}.";
            ls -l ${DBNSFP_DIR}/dbNSFP.gz*;
        fi
        echo;
        echo "Reinstalling dbNSFP may take a day."
        read -p  "Are you sure you want to download and reinstall dbNSFP? y [n] " -n 1 -r
        if [[ $REPLY =~ ^[Yy]$ ]]
        then
            log "downloading and installing dbNSFP";
            bash ${HERE}/download_dbnsfp_data.sh ${DOWNLOAD_DIR}/dbnsfp ${DOWNLOAD_DBNSFP_VERSION} ${DBNSFP_DIR}
        fi
    fi
    
    #FATHMM plugin
    if [[ DOWNLOAD_FATHMM == 1 ]]; then
        #update only needed if FATHMM changes
        FATHMM_DIR=${DOWNLOAD_DIR}/fathmm;
        mkdir - p ${FATHMM_DIR};
        log "downloading FATHMM data";
        bash ${HERE}/download_fathmm_data.sh ${FATHMM_DIR};
    
        log "creating fathmm database"
        #load fathmm database (only if fathmm was updated)
        ${MYSQL} -e "DROP DATABASE IF EXISTS fathmm";
        ${MYSQL} -e "CREATE SCHEMA IF NOT EXISTS fathmm";
        
        log "importing FATHMM database"
        ${MYSQL} -e "GRANT SELECT ON fathmm."'*'" TO '${VAX_USER}'@'%' IDENTIFIED BY '${VAX_PW}'";
        gunzip < ${FATHMM_DIR}/fathmm.v*.SQL.gz | ${MYSQL} -Dfathmm;
    
        log "creating ${PLUGINS_DIR}/fathmm"
        mkdir -p ${PLUGINS_DIR}/fathmm;
        log "copying ${FATHMM_DIR}/fathmm.py script to ${PLUGINS_DIR}/fathmm";
        cp -a ${FATHMM_DIR}/fathmm.py ${PLUGINS_DIR}/fathmm;
    
        #creating config file for downloaded FATHMM plugin (not part of VAX)
        log "creating ${PLUGINS_DIR}/fathmm/config.ini"
        cat > ${PLUGINS_DIR}/fathmm/config.ini << EOF
[DATABASE]
HOST = ${HOST}
PORT = ${PORT}
USER = ${VAX_USER}
PASSWD = ${VAX_PW}
DB = fathmm
EOF
    fi
fi

################################################################################
#                                   UniProt                                    #
################################################################################

echo;
UNIPROT_DIR=${DOWNLOAD_DIR}/uniprot;
mkdir -p ${UNIPROT_DIR};

if [[ ${DOWNLOAD} == 1 ]]; then
    #download UniProt data
    log "downloading UniProt data";
    bash ${HERE}/download_uniprot_data.sh ${UNIPROT_DIR};
    
    #create UniProt tab-delimited files
    #this function requires a reference to SWISS::Knife in your PERL5LIB
    #http://sourceforge.net/projects/swissknife/files/latest/download
    #export PERL5LIB=${PERL5LIB}:/path/to/Swissknife_1.69/lib
    log "creating UniProt tab-delimited files; ignore any defined(@array) is deprecated warning";
    perl ${HERE}/uniprot2db.pl -dir ${UNIPROT_DIR};
    
    #create Uniprot table definitions
    log "creating uniprot_human_protein_feature table definition";
    python ${HERE}/mysql_table.py -i ${UNIPROT_DIR}/uniprot_human_protein_feature.txt -t uniprot_human_protein_feature;
    log "creating uniprot_human_protein table definition";
    python ${HERE}/mysql_table.py -i ${UNIPROT_DIR}/uniprot_human_protein.txt -t uniprot_human_protein;
    log "creating uniprot_human_xref table definition";
    python ${HERE}/mysql_table.py -i ${UNIPROT_DIR}/uniprot_human_xref.txt -t uniprot_human_xref;
fi

if [[ ${INSTALL_DB} == 1 ]]; then
    #create UniProt tables
    log "creating uniprot_human_protein_feature mysql table";
    ${MYSQL} ${DATABASE} < ${UNIPROT_DIR}/uniprot_human_protein_feature.mysql
    log "creating uniprot_human_protein mysql table";
    ${MYSQL} ${DATABASE} < ${UNIPROT_DIR}/uniprot_human_protein.mysql
    log "creating uniprot_human_xref mysql table";
    ${MYSQL} ${DATABASE} < ${UNIPROT_DIR}/uniprot_human_xref.mysql
    
    #load UniProt data
    log "importing uniprot_human_protein_feature table";
    ${MYSQLIMPORT} --delete --local --ignore-lines=1 ${DATABASE} ${UNIPROT_DIR}/uniprot_human_protein_feature.txt;
    log "importing uniprot_human_protein table";
    ${MYSQLIMPORT} --delete --local --ignore-lines=1 ${DATABASE} ${UNIPROT_DIR}/uniprot_human_protein.txt;
    log "importing uniprot_human_xref table";
    ${MYSQLIMPORT} --delete --local --ignore-lines=1 ${DATABASE} ${UNIPROT_DIR}/uniprot_human_xref.txt;
    
    #index UniProt tables
    log "indexing uniprot_human_xref table";
    ${MYSQL} -e "ALTER TABLE ${DATABASE}.uniprot_human_xref
    ADD INDEX IX_RESOURCE_IDENTIFIER (RESOURCE_IDENTIFIER ASC),
    ADD INDEX IX_RESOURCE_ABBREVIATION (RESOURCE_ABBREVIATION ASC),
    ADD INDEX IX_UniProtKB_AC (UniProtKB_AC ASC);";
    
    #create enst_uniprot tables for vax access
    log "creating enst_uniprot table";
    ${MYSQL} -e "USE ${DATABASE};
    DROP TABLE IF EXISTS enst_uniprot;
    CREATE TABLE enst_uniprot
    (KEY IX_ENST_UniProtKB_AC_topic(ENST,UniProtKB_AC,topic),
    KEY IX_ENST_topic (ENST,topic),
    KEY IX_UniProtKB_AC_topic (UniProtKB_AC,topic))
    ENGINE=MyISAM DEFAULT CHARSET=latin1
    SELECT CAST(xe.RESOURCE_IDENTIFIER AS CHAR(15)) ENST, u.*
    FROM uniprot_human_protein u
    JOIN uniprot_human_xref xe
    ON u.UniProtKB_AC = xe.UniProtKB_AC AND xe.RESOURCE_ABBREVIATION = 'Ensembl' AND xe.RESOURCE_IDENTIFIER LIKE 'ENST%';";
    
    log "creating enst_uniprot_feature table";
    ${MYSQL} -e "USE ${DATABASE};
    DROP TABLE IF EXISTS enst_uniprot_feature;
    CREATE TABLE enst_uniprot_feature
    (KEY IX_ENST_feature_aaStart_aaEnd (ENST,feature,aaStart,aaEnd),
    KEY IX_UniProtKB_AC_feature_aaStart_aaEnd (UniProtKB_AC,feature,aaStart,aaEnd))
    ENGINE=MyISAM DEFAULT CHARSET=latin1
    SELECT DISTINCT CAST(xe.RESOURCE_IDENTIFIER AS CHAR(15)) ENST, u.*
    FROM uniprot_human_protein_feature u
    JOIN uniprot_human_xref xe
    ON u.UniProtKB_AC = xe.UniProtKB_AC AND xe.RESOURCE_ABBREVIATION = 'Ensembl' AND xe.RESOURCE_IDENTIFIER LIKE 'ENST%';";
    
    #create UniProt stored procedures for vax access
    log "creating enst2uniprot stored procedure";
    ${MYSQL} -e "USE ${DATABASE};
    DROP PROCEDURE IF EXISTS enst2uniprot;
    DELIMITER \$\$
    CREATE DEFINER=CURRENT_USER PROCEDURE enst2uniprot(enst CHAR(15))
    BEGIN
    SELECT u.topic, u.value
    FROM ${DATABASE}.enst_uniprot AS u
    WHERE u.ENST = enst
    AND COALESCE(u.value, '') <> ''
    ORDER BY u.topic;
    END\$\$
    DELIMITER ;";
    
    log "creating enst2uniprot_feature stored procedure";
    ${MYSQL} -e "USE ${DATABASE};
    DELIMITER \$\$
    DROP PROCEDURE IF EXISTS enst2uniprot_feature;
    CREATE DEFINER=CURRENT_USER PROCEDURE enst2uniprot_feature(enst CHAR(15), aaStart INT, aaEnd INT)
    BEGIN
    SELECT u.feature, u.aaStart, u.aaEnd, u.description
    FROM ${DATABASE}.enst_uniprot_feature AS u
    WHERE u.ENST = enst
    AND ((aaStart BETWEEN u.aaStart and u.aaEnd)
    OR  (aaEnd BETWEEN u.aaStart and u.aaEnd))
    ORDER BY u.feature, u.aaStart, u.aaEnd;
    END\$\$
    DELIMITER ;";
fi

################################################################################
#                                     KEGG                                     #
################################################################################

#KEGG installation requires the previous UniProt installation

echo;
KEGG_DIR=${DOWNLOAD_DIR}/kegg;
mkdir -p ${KEGG_DIR};

#download KEGG data
if [[ ${DOWNLOAD} == 1 ]]; then
    log "downloading KEGG data";
    perl ${HERE}/download_kegg_data.pl -o ${KEGG_DIR}/kegg_gene_pathways.txt;
    
    #create KEGG table definition
    log "creating kegg_gene_pathways table definition";
    python ${HERE}/mysql_table.py -i ${KEGG_DIR}/kegg_gene_pathways.txt -t kegg_gene_pathways;
fi

if [[ ${INSTALL_DB} == 1 ]]; then
    #create KEGG table
    log "creating kegg_gene_pathways mysql table";
    ${MYSQL} ${DATABASE} < ${KEGG_DIR}/kegg_gene_pathways.mysql
    
    #load KEGG data
    log "importing kegg_gene_pathways table";
    ${MYSQLIMPORT} --delete --local --ignore-lines=1 ${DATABASE} ${KEGG_DIR}/kegg_gene_pathways.txt;
    
    #index KEGG tables
    log "indexing kegg_gene_pathways table";
    ${MYSQL} -e "ALTER TABLE ${DATABASE}.kegg_gene_pathways
    ADD INDEX IX_gene_id (gene_id ASC);";
    
    #create ensg_kegg table for vax access
    log "creating ensg_kegg table";
    ${MYSQL} -e "USE ${DATABASE};
    DROP TABLE IF EXISTS ensg_kegg;
    CREATE TABLE ensg_kegg
    (PRIMARY KEY (ENSG,gene_id,path_id))
    ENGINE=MyISAM DEFAULT CHARSET=latin1
    SELECT DISTINCT CAST(xe.OPTIONAL_INFORMATION_2 AS CHAR(15)) ENSG, k.gene_id, k.path_id, LTRIM(RTRIM(k.pathway)) pathway
    FROM kegg_gene_pathways AS k
    JOIN uniprot_human_xref AS x
    ON k.gene_id = x.RESOURCE_IDENTIFIER AND x.RESOURCE_ABBREVIATION = 'KEGG'
    JOIN uniprot_human_xref AS xe
    ON x.UniProtKB_AC = xe.UniProtKB_AC AND xe.RESOURCE_ABBREVIATION = 'Ensembl' AND xe.OPTIONAL_INFORMATION_2 LIKE 'ENSG%'
    WHERE IFNULL(LTRIM(RTRIM(k.pathway)), '') <> '';";
    
    #create KEGG stored procedure for vax access
    log "creating ensg2kegg stored procedure";
    ${MYSQL} -e "USE ${DATABASE};
    DROP PROCEDURE IF EXISTS ensg2kegg;
    DELIMITER \$\$
    CREATE DEFINER=CURRENT_USER PROCEDURE ensg2kegg(ensg CHAR(15))
    BEGIN
    SELECT DISTINCT k.pathway
    FROM ${DATABASE}.ensg_kegg AS k
    WHERE k.ENSG = ensg
    ORDER BY k.path_id;
    END\$\$
    DELIMITER ;";
fi

################################################################################
#                             Human Protein Atlas                              #
################################################################################

echo;
HPA_DIR=${DOWNLOAD_DIR}/hpa;
mkdir -p ${HPA_DIR};


if [[ ${DOWNLOAD} == 1 ]]; then
    #download HPA data
    log "downloading HPA data";
    bash ${HERE}/download_hpa_data.sh  ${HPA_DIR};
    
    #create HPA tab-delimited files
    log "creating hpa_normal_tissue tab-delimited file"
    python ${HERE}/csv2tab.py -i ${HPA_DIR}/normal_tissue.csv -o ${HPA_DIR}/hpa_normal_tissue.tab
    log "creating hpa_subcellular_location tab-delimited file"
    python ${HERE}/csv2tab.py -i ${HPA_DIR}/subcellular_location.csv -o ${HPA_DIR}/hpa_subcellular_location.tab
    
    #create HPA table definitions
    log "creating hpa_normal_tissue table definition";
    python ${HERE}/mysql_table.py -i ${HPA_DIR}/hpa_normal_tissue.tab -t hpa_normal_tissue --header_line 1;
    log "creating hpa_subcellular_location table definition";
    python ${HERE}/mysql_table.py -i ${HPA_DIR}/hpa_subcellular_location.tab -t hpa_subcellular_location --header_line 1;
fi

if [[ ${INSTALL_DB} == 1 ]]; then
    #create HPA tables
    log "creating hpa_normal_tissue mysql table";
    ${MYSQL} ${DATABASE} < ${HPA_DIR}/hpa_normal_tissue.mysql;
    log "creating hpa_subcellular_location mysql table";
    ${MYSQL} ${DATABASE} < ${HPA_DIR}/hpa_subcellular_location.mysql;
    
    #load HPA data
    log "importing hpa_normal_tissue table";
    ${MYSQLIMPORT} --delete --local --ignore-lines=1 ${DATABASE} ${HPA_DIR}/hpa_normal_tissue.tab;
    log "importing hpa_subcellular_location table";
    ${MYSQLIMPORT} --delete --local --ignore-lines=1 ${DATABASE} ${HPA_DIR}/hpa_subcellular_location.tab;
    
    #index HPA tables
    log "indexing hpa_normal_tissue table";
    ${MYSQL} -e "ALTER TABLE ${DATABASE}.hpa_normal_tissue
    CHANGE COLUMN Gene Gene VARCHAR(15) NOT NULL,
    CHANGE COLUMN Tissue Tissue VARCHAR(17) NOT NULL,
    CHANGE COLUMN \`Cell type\` Cell_type VARCHAR(28) NOT NULL,
    CHANGE COLUMN \`Expression type\` Expression_type VARCHAR(8) NULL,
    ADD PRIMARY KEY (Gene, Tissue, Cell_type);";
   
    log "indexing hpa_subcellular_location table";
    ${MYSQL} -e "ALTER TABLE ${DATABASE}.hpa_subcellular_location
    CHANGE COLUMN Gene Gene VARCHAR(15) NOT NULL,
    CHANGE COLUMN \`Main location\` Main_location VARCHAR(81) NULL,
    CHANGE COLUMN \`Other location\` Other_location VARCHAR(86) NULL,
    CHANGE COLUMN \`Expression type\` Expression_type VARCHAR(8) NULL,
    ADD PRIMARY KEY (Gene);";
        
    #create HPA stored procedures for vax access
    log "creating hpa_tissue_cell_type_list stored procedure";
    ${MYSQL} -e "USE ${DATABASE};
    DROP PROCEDURE IF EXISTS hpa_tissue_cell_type_list;
    DELIMITER \$\$
    CREATE DEFINER=CURRENT_USER PROCEDURE hpa_tissue_cell_type_list()
    BEGIN
    SELECT DISTINCT h.Tissue, h.Cell_type
    FROM ${DATABASE}.hpa_normal_tissue AS h
    ORDER BY h.Tissue, h.Cell_type;
    END\$\$
    DELIMITER ;";
    
    log "creating ensg2hpa_subcellular_location stored procedure";
    ${MYSQL} -e "USE ${DATABASE};
    DROP PROCEDURE IF EXISTS ensg2hpa_subcellular_location;
    DELIMITER \$\$
    CREATE DEFINER=CURRENT_USER PROCEDURE ensg2hpa_subcellular_location(ensg CHAR(15))
    BEGIN
    SELECT DISTINCT h.Main_location AS location
    FROM ${DATABASE}.hpa_subcellular_location AS h
    WHERE h.Gene = ensg
    UNION
    SELECT DISTINCT h.Other_location AS location
    FROM ${DATABASE}.hpa_subcellular_location AS h
    WHERE h.Gene = ensg
    ORDER BY location;
    END\$\$
    DELIMITER ;";

    echo "creating ensg2hpa_tissue stored procedure";
    ${MYSQL} -e "USE ${DATABASE};
    DROP PROCEDURE IF EXISTS ensg2hpa_tissue;
    DELIMITER \$\$
    CREATE DEFINER=CURRENT_USER PROCEDURE ensg2hpa_tissue(ensg CHAR(15))
    BEGIN
    SELECT DISTINCT h.Tissue, h.Cell_type, h.Level
    FROM ${DATABASE}.hpa_normal_tissue AS h
    WHERE h.Gene = ensg
    ORDER BY h.Tissue, h.Cell_type;
    END\$\$
    DELIMITER ;";
fi

################################################################################
#                                  MitoCarta                                   #
################################################################################

echo;
#MITO_DIR=${DOWNLOAD_DIR}/mito;
#mkdir -p ${MITO_DIR};

#MitoCarta data and table definition are in the VAX utilities directory

if [[ ${INSTALL_DB} == 1 ]]; then
    #create MitoCarta table
    log "creating HumanMitoCartaAll mysql table";
    ${MYSQL} ${DATABASE} < ${HERE}/HumanMitoCartaAll.sql;
    
    #load MitoCarta data
    log "importing HumanMitoCartaAll table";
    ${MYSQLIMPORT} --delete --local --ignore-lines=1 ${DATABASE} ${HERE}/HumanMitoCartaAll.txt;
      
    #create MitoCarta stored procedure for vax access
    log "creating get_mitocarta_gene stored procedure";
    ${MYSQL} -e "USE ${DATABASE};
    DROP PROCEDURE IF EXISTS get_mitocarta_gene;
    DELIMITER \$\$
    CREATE DEFINER=CURRENT_USER PROCEDURE get_mitocarta_gene(hgnc varchar(50))
    BEGIN
    SELECT HumanMitoCartaAll.SYM AS mito_gene
    FROM ${DATABASE}.HumanMitoCartaAll
    WHERE HumanMitoCartaAll.SYM = hgnc AND MITOCARTA_LIST = 1;
    END\$\$
    DELIMITER ;";
fi

################################################################################
#                              Mouse Phenotypes                                #
################################################################################

echo;
MGI_DIR=${DOWNLOAD_DIR}/mgi;
mkdir -p ${MGI_DIR};


if [[ ${DOWNLOAD} == 1 ]]; then
    #download MGI data
    log "downloading MGI data";
    bash ${HERE}/download_mgi_data.sh  ${MGI_DIR};
    
    #create MGI tab-delimited files
    log "creating HMD_HumanPhenotype.rpt.flat and VOC_MammalianPhenotype.rpt.clean tab-delimited files."
    python ${HERE}/MGI_mouse_phenotype_files.py -H ${MGI_DIR}/HMD_HumanPhenotype.rpt -V ${MGI_DIR}/VOC_MammalianPhenotype.rpt
    
    #create MGI table definitions
    log "creating HMD_HumanPhenotype table definition";
    python ${HERE}/mysql_table.py -i ${MGI_DIR}/HMD_HumanPhenotype.rpt.flat -t HMD_HumanPhenotype;
    log "creating VOC_MammalianPhenotype table definition";
    python ${HERE}/mysql_table.py -i ${MGI_DIR}/VOC_MammalianPhenotype.rpt.clean -t VOC_MammalianPhenotype;
fi

if [[ ${INSTALL_DB} == 1 ]]; then
    #create MGI tables
    log "creating HMD_HumanPhenotype mysql table";
    ${MYSQL} ${DATABASE} < ${MGI_DIR}/HMD_HumanPhenotype.mysql;
    log "creating VOC_MammalianPhenotype mysql table";
    ${MYSQL} ${DATABASE} < ${MGI_DIR}/VOC_MammalianPhenotype.mysql;
    
    #load MGI data
    log "importing HMD_HumanPhenotype table";
    ${MYSQLIMPORT} --delete --local --ignore-lines=1 ${DATABASE} ${MGI_DIR}/HMD_HumanPhenotype.rpt.flat;
    log "importing VOC_MammalianPhenotype table";
    ${MYSQLIMPORT} --delete --local --ignore-lines=1 ${DATABASE} ${MGI_DIR}/VOC_MammalianPhenotype.rpt.clean;
    
    #index MGI tables
    log "indexing HMD_HumanPhenotype table";
    ${MYSQL} -e "ALTER TABLE ${DATABASE}.HMD_HumanPhenotype
    ADD INDEX IX_Mammalian_Phenotype_ID (Mammalian_Phenotype_ID ASC),
    ADD INDEX IX_Human_Gene (Human_Gene ASC);";

    log "indexing VOC_MammalianPhenotype table";
    ${MYSQL} -e "ALTER TABLE ${DATABASE}.VOC_MammalianPhenotype
    ADD INDEX IX_Mammalian_Phenotype_ID (Mammalian_Phenotype_ID ASC);";
        
    #create MGI stored procedure for vax access
    log "creating hgnc2mouse_phenotype stored procedure";
    ${MYSQL} -e "USE ${DATABASE};
    DROP PROCEDURE IF EXISTS hgnc2mouse_phenotype;
    DELIMITER \$\$
    CREATE DEFINER=CURRENT_USER PROCEDURE hgnc2mouse_phenotype(hgnc VARCHAR(15))
    BEGIN
    SELECT DISTINCT v.Mammalian_Phenotype_Name
    from vw.VOC_MammalianPhenotype as v
    join vw.HMD_HumanPhenotype as h
    on v.Mammalian_Phenotype_ID = h.Mammalian_Phenotype_ID
    WHERE h.Human_Gene = hgnc 
    ORDER BY v.Mammalian_Phenotype_Name;
    END\$\$
    DELIMITER ;";
fi

################################################################################
#                                     OMIM                                     #
################################################################################

echo;
OMIM_DIR=${DOWNLOAD_DIR}/omim;
mkdir -p ${OMIM_DIR};

if [[ ${DOWNLOAD} == 1 ]]; then
    #download OMIM data
    log "downloading OMIM data";
    bash ${HERE}/download_omim_data.sh  ${OMIM_DIR};
    
    #create OMIM tab-delimited files
    log "creating OMIM tab-delimited files"
    python ${HERE}/omim2db.py -d ${OMIM_DIR};
    
    #create OMIM table definitions
    log "creating omim_genemap table definition";
    python ${HERE}/mysql_table.py -i ${OMIM_DIR}/omim_genemap.txt -t omim_genemap;
    
fi

if [[ ${INSTALL_DB} == 1 ]]; then
    #create OMIM table
    log "creating omim_genemap mysql table";
    ${MYSQL} ${DATABASE} < ${OMIM_DIR}/omim_genemap.mysql;
    
    #load OMIM data
    log "importing omim_genemap table";
    ${MYSQLIMPORT} --delete --local --ignore-lines=1 ${DATABASE} ${OMIM_DIR}/omim_genemap.txt;
    
    #index OMIM tables
    log "indexing omim_genemap table";
    ${MYSQL} -e "ALTER TABLE ${DATABASE}.omim_genemap
    ADD INDEX IX_gene (gene ASC);";
    
    #create OMIM table for vax access
    #it is possible to add ENSG identifiers to the omim_genemap table
    #but here we will use the HGNC identifiers in column omim_genemap.gene
    
    #create OMIM stored procedure for vax access
    #a different procedure would be used if an ENSG column were available
    log "creating hgnc2omim stored procedure";
    ${MYSQL} -e "USE ${DATABASE};
    DROP PROCEDURE IF EXISTS hgnc2omim;
    DELIMITER \$\$
    CREATE DEFINER=CURRENT_USER PROCEDURE hgnc2omim(hgnc VARCHAR(255))
    BEGIN
    SELECT distinct CAST(CONCAT(o.disorder, ' (', o.mim_locus, ')') AS CHAR(255)) AS disorder
    FROM ${DATABASE}.omim_genemap AS o
    WHERE o.gene = hgnc
    ORDER BY o.disorder;
    END\$\$
    DELIMITER ;";
fi

################################################################################
#                                   RefSeq                                     #
################################################################################

echo;
REFSEQ_DIR=${DOWNLOAD_DIR}/refseq;
mkdir -p ${REFSEQ_DIR};

if [[ ${DOWNLOAD} == 1 ]]; then
    #download RefSeq data
    log "downloading RefSeq data";
    bash ${HERE}/download_refseq_data.sh  ${REFSEQ_DIR};
    
    #create RefSeq tab-delimited files
    log "creating RefSeq tab-delimited file";
    python ${HERE}/refgene2db.py -r ${REFSEQ_DIR}/gene_RefSeqGene -i ${REFSEQ_DIR}/refseqgene*.genomic.gbff.gz -o ${REFSEQ_DIR}/refseq_summary.txt
    
    #create RefSeq table definitions
    log "creating refseq_summary table definition";
    python ${HERE}/mysql_table.py -i ${REFSEQ_DIR}/refseq_summary.txt -t refseq_summary;
    
fi

if [[ ${INSTALL_DB} == 1 ]]; then
    #create refseq_summary table
    log "creating RefSeq mysql table";
    ${MYSQL} ${DATABASE} < ${REFSEQ_DIR}/refseq_summary.mysql;
    
    #load RefSeq data
    log "importing refseq_summary table";
    ${MYSQLIMPORT} --delete --local --ignore-lines=1 ${DATABASE} ${REFSEQ_DIR}/refseq_summary.txt;
    
    #index RefSeq tables
    log "indexing refseq_summary table";
    ${MYSQL} -e "ALTER TABLE ${DATABASE}.refseq_summary
    ADD INDEX IX_GeneSymbol (GeneSymbol ASC);";
    
    #create RefSeq table for vax access
    #it is possible to add ENSG identifiers to the refseq_summary table
    #but here we will use the HGNC identifiers in column refseq_summary.GeneSymbol
    
    #create RefSeq stored procedure for vax access
    #a different procedure would be used if an ENSG column were available
    log "creating hgnc2refseq_summary stored procedure";
    ${MYSQL} -e "USE ${DATABASE};
    DROP PROCEDURE IF EXISTS hgnc2refseq_summary;
    DELIMITER \$\$
    CREATE DEFINER=CURRENT_USER PROCEDURE hgnc2refseq_summary(hgnc VARCHAR(255))
    BEGIN
    SELECT DISTINCT r.SUMMARY AS summary
    FROM ${DATABASE}.refseq_summary AS r
    WHERE r.GeneSymbol = hgnc
    ORDER BY r.SUMMARY;
    END\$\$
    DELIMITER ;";
fi

echo;
log "VAX install done.";
echo;
echo "The downloaded data files in ${DOWNLOAD_DIR} may be deleted or saved for a future update."
echo;

exit;
