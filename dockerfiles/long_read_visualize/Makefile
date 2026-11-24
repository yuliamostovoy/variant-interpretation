VERSION = 0.3
#TAG1 = quay.io/ymostovoy/lr-visualize:$(VERSION)
#TAG2 = quay.io/ymostovoy/lr-visualize:test

TAG1 = us-central1-docker.pkg.dev/talkowski-sv-gnomad/ym-docker/lr-visualize:$(VERSION)
TAG2 = us-central1-docker.pkg.dev/talkowski-sv-gnomad/ym-docker/lr-visualize:test

all: build push

build:
	docker build -t $(TAG1) -t $(TAG2) .

push:
	docker push $(TAG1)
	docker push $(TAG2)
