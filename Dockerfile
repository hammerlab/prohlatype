FROM ocaml/opam:ubuntu-16.04_ocaml-4.05.0
RUN opam repo add main https://opam.ocaml.org
RUN opam update --yes
RUN sudo apt-get update --yes
RUN sudo apt-get install --yes git libgmp-dev zlib1g-dev
RUN git clone https://github.com/hammerlab/prohlatype.git
WORKDIR ./prohlatype
RUN opam info ppx_deriving
RUN opam info ppx_deriving_yojson
RUN opam pin add -k git prohlatype .
RUN opam config exec -- make build tools

