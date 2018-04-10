FROM ocaml/opam:ubuntu-16.04_ocaml-4.05.0
ARG prohlatype_tag
RUN opam repo add main https://opam.ocaml.org
RUN opam update --yes
RUN sudo apt-get update --yes
RUN sudo apt-get install --yes git libgmp-dev zlib1g-dev libffi-dev
RUN git clone https://github.com/hammerlab/prohlatype.git
WORKDIR ./prohlatype
RUN git checkout tags/$prohlatype_tag
RUN opam pin add -k git prohlatype .
RUN opam config exec -- make
RUN ln -s _build/default/src/app/align2fasta.exe align2fasta
RUN ln -s _build/default/src/app/mhc2gpdf.exe mhc2gpdf
RUN ln -s _build/default/src/app/mpjson2tsv.exe mpjson2tsv
RUN ln -s _build/default/src/app/multi_par.exe multi_par
RUN ln -s _build/default/src/app/par_type.exe par_type
