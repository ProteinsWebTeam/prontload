# pronto-loader

Refresh Pronto with the latest data from InterPro, GOA, and UniProt.

## Installation

Get the latest version of the code and install requirements (Python 3.4+):

```sh
$ git clone https://github.com/ProteinsWebTeam/pronto-loader.git
$ cd pronto-loader/
$ python setup.py install
```

## Configuration

Create or edit the `config.json` file to set the following settings:

* `dsn`: Oracle DSN in the `user/password@host:port/service` format.
* `schema`: Oracle schema where to create/refresh tables. Be sure that the user defined in `dsn` has the right priviles (e.g. `CREATE TABLE`, `DROP ANY TABLE`, etc.).
* `max_gap`: Used to group proteins by their match structure. Historically set to `20`.

## Run

The only requirement argument is the path the `config.json`.

```sh
pronto-update config.json [OPTIONS]
```

Available options are:

| Option        | Description                                                                    |
| ------------- |--------------------------------------------------------------------------------|
| -s, --step    | space-separated list of steps to run (default: all)                            |
| -p, --threads | number of child processes to use for the `matches` step (default: `2`)         |
| -t, --temp    | directory for temporary files (default: probably `/tmp/`)                      |
| -o, --output  | output file for the Swiss-Prot descriptions report for curators                |

### Steps

| Name          | Description                                                                    |
| ------------- |--------------------------------------------------------------------------------|
| databases     | Copy database info (e.g. name, version) from InterPro                          |
| synonyms      | Create synonyms for InterPro tables that are updated daily by curators         |
| comments      | Copy Swiss-Prot comments from UniProt                                          |
| descriptions  | Copy protein descriptions from UniProt                                         |
| enzymes       | Copy EC numbers ([ENZYME](https://enzyme.expasy.org/)) FROM UniProt            |
| protein2go    | Copy protein GO annotations from GOA                                           |
| methods       | Copy member database signatures from InterPro                                  |
| taxonomies    | Copy taxonomy data from InterPro                                               |
| terms         | Copy GO terms from GOA                                                         |
| proteins      | Copy proteins from InterPro                                                    |
| matches       | Copy protein matches from InterPro, and run predictions                        |
| report        | Export gained/lost protein descriptions per InterPro entries                   |
