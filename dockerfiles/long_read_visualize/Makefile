VERSION = 0.1
TAG1 = quay.io/ymostovoy/lr-visualize:$(VERSION)
TAG2 = quay.io/ymostovoy/lr-visualize:latest

all: build push

build:
	docker build -t $(TAG1) -t $(TAG2) .

push:
	docker push $(TAG1)
	docker push $(TAG2)
