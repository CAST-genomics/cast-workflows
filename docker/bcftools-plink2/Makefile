# Definitions
repository = gcr.io/ucsd-medicine-cast
identifier = bcftools-plink2
version = 0.0.0
git_commit ?= $(shell git log --pretty=oneline -n 1 | cut -f1 -d " ")
name = ${repository}/${identifier}
tag = ${version}--${git_commit}

# Steps
build:
	# do the docker build
	docker build --platform="linux/amd64" -t ${name}:${tag} .
	# also build "native" (mostly for new macs...)
	docker build -t ${name}:${tag}-native .
	docker tag ${name}:${tag} ${name}:latest
	docker tag ${name}:${tag}-native ${name}:latest-native

push:
	# Requires ~/.dockercfg
	docker push ${name}:${tag}
	docker push ${name}:latest-native
