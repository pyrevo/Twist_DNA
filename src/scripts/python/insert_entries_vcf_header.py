#!/bin/python3


with open(snakemake.output[0], 'w') as output:
    with open(snakemake.input[0]) as input:
        output.write(next(input))

        for entry in snakemake.params.entries:
            output.write("##" + snakemake.params.type + "=<" +
                         ",".join([e[0] + "=" + e[1] for e in entry]) + ">\n"
                         )

        for row in input:
            output.write(row)
