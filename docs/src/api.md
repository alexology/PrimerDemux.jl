# API

```@meta
CurrentModule = PrimerDemux
```

```@docs
PrimerDemux
```

## Processing files

```@docs
remove_primers
demultiplex
TrimStats
```

## Describing adapters and primers

```@docs
front_adapter
back_adapter
adapters_from_fasta
Adapter
SingleAdapter
LinkedAdapter
Placement
FRONT
BACK
PREFIX
SUFFIX
FRONT_NOT_INTERNAL
BACK_NOT_INTERNAL
```

## Single sequences

```@docs
find_adapter
AdapterMatch
trim_sequence
reverse_complement
quality_trim_index
```
