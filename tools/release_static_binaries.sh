# in
# docker run --rm -ti ocaml/opam:alpine_ocaml-4.06.0 bash

export PTAG="0.9.1"

cd opam-repository/
git pull origin master

export OPAMYES="true"
sudo apk add m4 pkgconfig libffi-dev zlib-dev
opam update && opam upgrade
opam depext ssl
cd ..
git clone https://github.com/hammerlab/prohlatype.git
cd prohlatype
git checkout tags/${PTAG}
opam pin add -k git prohlatype .
make release
mkdir zip_me
cp _build/default/src/app/*.exe zip_me/
cd zip_me
zip prohlatype.${PTAG}.zip *.exe

# outside:
docker cp ${CONTAINER}:/home/opam/prohlatype/zip_me/prohlatype.${PTAG}.zip .
