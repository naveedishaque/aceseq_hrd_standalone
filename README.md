# ACEseq Workflow HRD calculator plugin

This workflow runs the calculation of the HRD step for ACEseq as a standalone plugin. This is based on the ACEseq (TODO: add link to hash ref on github)

## Pre-requisites

### conda environment

Setup a conda environment: `conda create --name hrd python=2.7 numpy bedtools r=3.3`

On the eils-hpc this env can be sourced here: `/home/ishaquen/miniconda3/bin/activate /home/ishaquen/miniconda3/envs/hrd`. **This is the default action for loading the env!**


## Usage

```
  sh estimateHRDScore.sh [/path/to/ACEseq_roddy_runtimeconfig.sh] [/path/to/analysis/ACEseq_dir] [pid]

```


## Versioning

We use [SemVer](http://semver.org/) for versioning. For the versions available, see the [tags on this repository](https://github.com/your/project/tags).
 - v1.0.0: first working version

## Authors

Naveed Ishaque

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details
