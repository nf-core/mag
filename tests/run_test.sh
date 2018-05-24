#!/usr/bin/env bash

# print_usage()
function print_usage {
  echo -e  "\nUsage:\t$0\n" \
    "\t\t[-b (example flag)\n" \
    "\t\t[-t <test data directory>]\n" \
    "\t\t[-d <docker image>]\n" \
    "\t\t[-s <singularity image>]\n" \
    "\t\t[-h (show this help message)]" >&2 ;
}

# Check that we have required commands
curl --version >/dev/null 2>&1 || { echo >&2 "I require curl, but it's not installed. Aborting."; exit 1; }
tar --version >/dev/null 2>&1 || { echo >&2 "I require tar, but it's not installed. Aborting."; exit 1; }
nextflow -v >/dev/null 2>&1 || { echo >&2 "I require nextflow, but it's not installed. If you hava Java, run 'curl -fsSL get.nextflow.io | bash'. If not, install Java."; exit 1; }

# Detect Travis fork for dockerhub image if we can
dockerfl=""
if [[ ! -z "$TRAVIS_REPO_SLUG" ]]; then
    dockerimg=$(echo "$TRAVIS_REPO_SLUG" | awk '{print tolower($0)}')
    echo "Detected repo as '$TRAVIS_REPO_SLUG' - using docker image '$dockerimg'"
    dockerfl="-with-docker $dockerimg"
fi

# Look for an existing test data directory
data_path="/tmp"
if [ -d "./test_data" ]
then
    data_path="./test_data"
    echo "Found data directory in current working directory, using ./test_data/"
fi
data_dir=${data_path}/nf-core/mag_test_set

# command line options
pipelinescript="../main.nf"
profile="-profile docker --max_cpus 2 --max_memory '7.GB' --max_time '48.h'"
singularityfl=""

while getopts ":brnpuht:d:s:" opt; do
  case $opt in
    b)
      echo "This is an example flag, it doesn't do anything" >&2
      ;;
    t)
      echo "Test data path specified" >&2
      data_path=$OPTARG
      ;;
    d)
      echo "Using docker image $OPTARG" >&2
      dockerfl="-with-docker $OPTARG"
      ;;
    s)
      echo "Using singularity image $OPTARG" >&2
      singularityfl="-with-singularity $OPTARG"
      ;;
    h)
      print_usage
      exit
      ;;
    :)
      echo -e "\nOption -$OPTARG requires an argument." >&2
      print_usage
      exit 1;
      ;;
    \?)
      print_usage
      echo -e "\nInvalid option: -$OPTARG" >&2
      exit 1
      ;;
  esac
done


if [ -d $data_dir ]
then
    echo "Found existing test set, using $data_dir"
else
    echo "Downloading test set..."
    curl -J -L https://github.com/HadrienG/test-datasets/raw/mag/test_data/mag_test_set.tar.bz2 > ${data_path}/mag_test_set.tar.bz2
    echo "Unpacking test set..."
    tar xvjf ${data_path}/mag_test_set.tar.bz2 -C ${data_path}
    echo "Done"
fi

# Run name
run_name="Test nf-core/mag Run: "$(date +%s)
# -name \"$run_name\"
cmd="nextflow run $pipelinescript -resume $profile $dockerfl $singularityfl --reads \"${data_path}/*R{1,2}.fastq.gz\""
echo "Starting nextflow... Command:"
echo $cmd
echo "-------------------------------------------------------"
eval $cmd
