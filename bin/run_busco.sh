#! /usr/bin/env bash

# Originally written by Sabrina Krakau and James Fellows Yates and released
# under the MIT license.
# See git repository (https://github.com/nf-core/mag) for full license text.

p=$1
cp_augustus_config=$2
db=$3
bin=$4
task_cpus=$5
lineage_dataset_provided=$6
busco_clean=$7
extra_args=$8

# ensure augustus has write access to config directory
if [ ${cp_augustus_config} = "Y" ]; then
    cp -r /usr/local/config/ augustus_config/
    export AUGUSTUS_CONFIG_PATH=augustus_config
fi

# place db in extra folder to ensure BUSCO recognizes it as path (instead of downloading it)
if [ ${lineage_dataset_provided} = "Y" ]; then
    mkdir dataset
    mv ${db} dataset/
fi

# set nullgob: if pattern matches no files, expand to a null string rather than to itself
shopt -s nullglob

# only used for saving busco downloads
most_spec_db="NA"

if
    busco ${p} \
        --mode genome \
        --in ${bin} \
        --cpu ${task_cpus} \
        ${extra_args} \
    --out "BUSCO" >${bin}_busco.log 2>${bin}_busco.err
then

    # get name of used specific lineage dataset
    summaries=(BUSCO/short_summary.specific.*.BUSCO.txt)
    if [ ${#summaries[@]} -ne 1 ]; then
        echo "ERROR: none or multiple 'BUSCO/short_summary.specific.*.BUSCO.txt' files found. Expected one."
        exit 1
    fi
    [[ $summaries =~ BUSCO/short_summary.specific.(.*).BUSCO.txt ]]
    db_name_spec="${BASH_REMATCH[1]}"
    most_spec_db=${db_name_spec}
    echo "Used specific lineage dataset: ${db_name_spec}"

    if [ ${lineage_dataset_provided} = "Y" ]; then
        cp BUSCO/short_summary.specific.${db_name_spec}.BUSCO.txt short_summary.specific_lineage.${db_name_spec}.${bin}.txt

        # if lineage dataset is provided, BUSCO analysis does not fail in case no genes can be found as when using the auto selection setting
        # report bin as failed to allow consistent warnings within the pipeline for both settings
        if egrep -q $'WARNING:\tBUSCO did not find any match.' ${bin}_busco.log; then
            echo "WARNING: BUSCO could not find any genes for the provided lineage dataset! See also ${bin}_busco.log."
            echo -e "${bin}\tNo genes" >"${bin}_busco.failed_bin.txt"
        fi
    else
        # auto lineage selection
        if { egrep -q $'INFO:\t\\S+ selected' ${bin}_busco.log &&
            egrep -q $'INFO:\tLineage \\S+ is selected, supported by ' ${bin}_busco.log; } ||
            { egrep -q $'INFO:\t\\S+ selected' ${bin}_busco.log &&
                egrep -q $'INFO:\tThe results from the Prodigal gene predictor indicate that your data belongs to the mollicutes clade. Testing subclades...' ${bin}_busco.log &&
                egrep -q $'INFO:\tUsing local lineages directory ' ${bin}_busco.log; }; then
            # the second statement is necessary, because certain mollicute clades use a different genetic code, are not part of the BUSCO placement tree, are tested separately
            # and cause different log messages
            echo "Domain and specific lineage could be selected by BUSCO."
            cp BUSCO/short_summary.specific.${db_name_spec}.BUSCO.txt short_summary.specific_lineage.${db_name_spec}.${bin}.txt

            db_name_gen=""
            summaries_gen=(BUSCO/short_summary.generic.*.BUSCO.txt)
            if [ ${#summaries_gen[@]} -lt 1 ]; then
                echo "No 'BUSCO/short_summary.generic.*.BUSCO.txt' file found. Assuming selected domain and specific lineages are the same."
                cp BUSCO/short_summary.specific.${db_name_spec}.BUSCO.txt short_summary.domain.${db_name_spec}.${bin}.txt
                db_name_gen=${db_name_spec}
            else
                [[ $summaries_gen =~ BUSCO/short_summary.generic.(.*).BUSCO.txt ]]
                db_name_gen="${BASH_REMATCH[1]}"
                echo "Used generic lineage dataset: ${db_name_gen}"
                cp BUSCO/short_summary.generic.${db_name_gen}.BUSCO.txt short_summary.domain.${db_name_gen}.${bin}.txt
            fi

            for f in BUSCO/run_${db_name_gen}/busco_sequences/single_copy_busco_sequences/*faa; do
                cat BUSCO/run_${db_name_gen}/busco_sequences/single_copy_busco_sequences/*faa | gzip >${bin}_buscos.${db_name_gen}.faa.gz
                break
            done
            for f in BUSCO/run_${db_name_gen}/busco_sequences/single_copy_busco_sequences/*fna; do
                cat BUSCO/run_${db_name_gen}/busco_sequences/single_copy_busco_sequences/*fna | gzip >${bin}_buscos.${db_name_gen}.fna.gz
                break
            done

        elif egrep -q $'INFO:\t\\S+ selected' ${bin}_busco.log && egrep -q $'INFO:\tNo marker genes were found. Root lineage \\S+ is kept' ${bin}_busco.log; then
            echo "Domain could be selected by BUSCO, but no more specific lineage."
            cp BUSCO/short_summary.specific.${db_name_spec}.BUSCO.txt short_summary.domain.${db_name_spec}.${bin}.txt

        elif egrep -q $'INFO:\t\\S+ selected' ${bin}_busco.log && egrep -q $'INFO:\tNot enough markers were placed on the tree \\([0-9]*\\). Root lineage \\S+ is kept' ${bin}_busco.log; then
            echo "Domain could be selected by BUSCO, but no more specific lineage."
            cp BUSCO/short_summary.specific.${db_name_spec}.BUSCO.txt short_summary.domain.${db_name_spec}.${bin}.txt

        elif egrep -q $'INFO:\t\\S+ selected' ${bin}_busco.log && egrep -q $'INFO:\tRunning virus detection pipeline' ${bin}_busco.log; then
            # TODO double-check if selected dataset is not one of bacteria_*, archaea_*, eukaryota_*?
            echo "Domain could not be selected by BUSCO, but virus dataset was selected."
            cp BUSCO/short_summary.specific.${db_name_spec}.BUSCO.txt short_summary.specific_lineage.${db_name_spec}.${bin}.txt
        else
            echo "ERROR: Some not expected case occurred! See ${bin}_busco.log." >&2
            exit 1
        fi
    fi

    for f in BUSCO/run_${db_name_spec}/busco_sequences/single_copy_busco_sequences/*faa; do
        cat BUSCO/run_${db_name_spec}/busco_sequences/single_copy_busco_sequences/*faa | gzip >${bin}_buscos.${db_name_spec}.faa.gz
        break
    done
    for f in BUSCO/run_${db_name_spec}/busco_sequences/single_copy_busco_sequences/*fna; do
        cat BUSCO/run_${db_name_spec}/busco_sequences/single_copy_busco_sequences/*fna | gzip >${bin}_buscos.${db_name_spec}.fna.gz
        break
    done

elif egrep -q $'ERROR:\tNo genes were recognized by BUSCO' ${bin}_busco.err; then
    echo "WARNING: BUSCO analysis failed due to no recognized genes! See also ${bin}_busco.err."
    echo -e "${bin}\tNo genes" >"${bin}_busco.failed_bin.txt"

elif egrep -q $'INFO:\t\\S+ selected' ${bin}_busco.log && egrep -q $'ERROR:\tPlacements failed' ${bin}_busco.err; then
    echo "WARNING: BUSCO analysis failed due to failed placements! See also ${bin}_busco.err. Still using results for selected generic lineage dataset."
    echo -e "${bin}\tPlacements failed" >"${bin}_busco.failed_bin.txt"

    message=$(egrep $'INFO:\t\\S+ selected' ${bin}_busco.log)
    [[ $message =~ INFO:[[:space:]]([_[:alnum:]]+)[[:space:]]selected ]]
    db_name_gen="${BASH_REMATCH[1]}"
    most_spec_db=${db_name_gen}
    echo "Used generic lineage dataset: ${db_name_gen}"
    cp BUSCO/auto_lineage/run_${db_name_gen}/short_summary.txt short_summary.domain.${db_name_gen}.${bin}.txt

    for f in BUSCO/auto_lineage/run_${db_name_gen}/busco_sequences/single_copy_busco_sequences/*faa; do
        cat BUSCO/auto_lineage/run_${db_name_gen}/busco_sequences/single_copy_busco_sequences/*faa | gzip >${bin}_buscos.${db_name_gen}.faa.gz
        break
    done
    for f in BUSCO/auto_lineage/run_${db_name_gen}/busco_sequences/single_copy_busco_sequences/*fna; do
        cat BUSCO/auto_lineage/run_${db_name_gen}/busco_sequences/single_copy_busco_sequences/*fna | gzip >${bin}_buscos.${db_name_gen}.fna.gz
        break
    done

else
    echo "ERROR: BUSCO analysis failed for some unknown reason! See also ${bin}_busco.err." >&2
    exit 1
fi

# additionally output genes predicted with Prodigal (GFF3)
if [ -f BUSCO/logs/prodigal_out.log ]; then
    mv BUSCO/logs/prodigal_out.log "${bin}_prodigal.gff"
fi

# output value of most_spec_db
echo ${most_spec_db} >info_most_spec_db.txt

# if needed delete temporary BUSCO files
if [ ${busco_clean} = "Y" ]; then
    find . -depth -type d -name "augustus_config" -execdir rm -rf "{}" \;
    find . -depth -type d -name "auto_lineage" -execdir rm -rf "{}" \;
    find . -depth -type d -name "run_*" -execdir rm -rf "{}" +
fi
