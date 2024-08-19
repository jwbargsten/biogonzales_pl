.PHONY: all

ifndef ANALYSIS_VERSION
  ifneq ("$(wildcard av)","")
    ANALYSIS_VERSION=$(shell cat av)
  else
    ANALYSIS_VERSION="."
  endif
endif

$(info analysis version: $(ANALYSIS_VERSION))
export ANALYSIS_VERSION

AV=$(ANALYSIS_VERSION)
NAME=$(lastword $(subst /, ,$(abspath .)))

all: dir

dir:
	@[ -d $(AV) ] || mkdir $(AV)

help:
	@grep -E '^[a-z0-9A-Z_.-]+:.*?## .*$$' $(MAKEFILE_LIST) | awk 'BEGIN {FS = ":.*?## "}; {printf "\033[1m%-30s\033[0m %s\n", $$1, $$2}'

