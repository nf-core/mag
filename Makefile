clean:
	rm -rf work/
	rm -rf .nextflow.log*
	rm -rf tests/

lint:
	nf-core lint .

test: clean lint
	mkdir -p tests/data && cd tests/data &&\
	curl -O -J -L https://github.com/HadrienG/test-datasets/raw/mag/test_data/test_mag_R1.fastq.gz &&\
	curl -O -J -L https://github.com/HadrienG/test-datasets/raw/mag/test_data/test_mag_R2.fastq.gz &&\
	cd .. && nextflow run ../main.nf -profile test,docker

docker:
	docker rmi --force hadrieng/mag:0.1.0dev
	docker build -t hadrieng/mag:0.1.0dev .
