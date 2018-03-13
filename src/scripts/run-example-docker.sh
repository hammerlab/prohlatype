#! /bin/sh

help () {
    cat <<'EOF'

Usage: sh src/scripts/run-example-docker.sh [prepare]

This script performs an example end-to-end run as described in the
`README.md` file.

It should require only `docker` (configured to run without `sudo`),
`wget`, and `git`; it fetches the data and runs everything in a
docker container.

The data is downloaded an computed in the `$PWD/tmp/` directory.

Bonus: sh src/scripts/run-example-docker.sh interactive
   drops you in a bash shell with all Prohlatype tools installed, as
   well as `bwa` and `samtools` (1.7).

EOF
}

set -e

mothership_docker_image=leonidr/prohlatype:0.8.0
docker_image=local/prohlatest:latest
tmp=$PWD/tmp/
sample_uri=https://raw.githubusercontent.com/hammerlab/prohlatype/master/tools/test_reads.fastq
prepare () {
    docker pull $mothership_docker_image

    mkdir -p $tmp
    cat > $tmp/Dockerfile <<EOF
FROM $mothership_docker_image
RUN sudo apt-get install -y bwa wget libarchive-dev libbz2-dev liblzma-dev
RUN wget https://github.com/samtools/samtools/releases/download/1.7/samtools-1.7.tar.bz2
RUN tar xvfj samtools-1.7.tar.bz2
RUN bash -c "cd samtools-1.7 && ./configure && make && sudo make install"
EOF

    docker build -t $docker_image $tmp/


    if [ -f imgthla/hla.dat ] ; then
        echo "IMGT HLA already there"
    else
        git clone https://github.com/ANHIG/IMGTHLA.git imgthla
    fi
    rm -fr tmp/results
    mkdir tmp/results
    if [ -f tmp/results/sample.fastq ] ; then
        echo "Test sample already there"
    else
        wget $sample_uri -O tmp/results/sample.fastq
    fi
    chmod -R 777 imgthla tmp
}

docker_options="-v $PWD/imgthla:/imgthla -v $PWD/tmp/results/:/results"
interactive () {
    docker run -it $docker_options $docker_image bash
}

runcmd () {
    local cmd="$1"
    echo "/\\__/o< $1"
    docker run $docker_options $docker_image bash -c "umask 000 ; $cmd"
}

default () {
    prepare
    runcmd "align2fasta --version"
    runcmd "multi_par --version"
    runcmd "bwa 2>&1 | grep Version"
    runcmd "samtools 2>&1 | grep Version"
    runcmd "find /imgthla -type f | wc -l"
    runcmd "align2fasta /imgthla/alignments -o /results/imputed_hla_class_I"
    runcmd "bwa index /results/imputed_hla_class_I.fasta"
    runcmd "bwa mem /results/imputed_hla_class_I.fasta /results/sample.fastq > /results/sample.sam"
    runcmd "samtools view -F 4 -b -T /results/imputed_hla_class_I.fasta /results/sample.sam -o /results/sample.bam"
    runcmd "samtools fastq /results/sample.bam > /results/sample_filtered.fastq"
    runcmd "multi_par /imgthla/alignments /results/sample_filtered.fastq -o /results/sample_output.tsv"
}

if [ "$1" = "" ] ; then
    default
else
    $*
fi
