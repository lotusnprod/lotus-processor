# Makefile README

This is a quick guide concerning the use of `make`.

## Instructions

### Mode

First of all, you can choose between the minimal working example `min` or making the whole DB `full`. By defaut make will run the `min` example for tests. To override it, simply do:

```console
make MODE=full <target>
```

### Docker

The `docker-build` and `docker-bash` commands are for our CI integration, you normally do not need them.

### Gathering the data

First, you need to perform the `gathering` of your initial data. To do so, you have multiple options:

- `gathering-full` will perform **all** gathering steps (**long**)
- `gathering-databases-full` will perform **all** databases related gathering steps (**long**). For a quicker version once you performed gathering all sub-steps simply do `gathering-databases`.
- `gathering-translation-full` will perform **all** translation dictionaries related gathering steps.

### Curating the data

Once done, you can do the `curating` steps of your initial data. If you are in a hurry, following command will do everything:

```console
make curating
```

Globally, `curating` is divided in three steps:

- `curating-1-integrating`: integrates all previously gathered data
- `curating-editing`: will perform the editing (cleaning) of the data
- `curating-3-integrating`: will perform the integration of cleaned data with initial data

Globally, `curating-editing` is also divided in three steps:

- `curating-editing-structure`: will perform the curation of the structures
- `curating-editing-organism`: will perform the curation of the organisms
- `curating-editing-reference`: will perform the curation of the references

Each curation step for each object is divided in more substeps. Fore more advanced details please see [the makefile](../Makefile).

### Analysing the data

After curation, results can be analysed:

```console
make analysing
```

### Visualizing the data

As a final step, figures can be reproduced using:

```console
make visualizing
```

### tl;dr

To perform all steps except gathering of the data in one shot:

```console
make curating-and-analysing-and-visualizing
```

additional details about the command arguments can be found [here](https://www.gnu.org/software/make/manual/make.html#Options-Summary)

## More info

If you want to know more simply visit <https://www.gnu.org/software/make/manual/make.html>