# Troubleshooting

## Running `domtblop.py parse <domtblout>` raises a ValueError from Biopython

If you see an error like this:

```
ValueError: The ID or alternative IDs of Hit '<hmm name>' exists in this QueryResult.
```

It is most likely to be caused by query sequences sharing the same name.
If they have the same sequence, one possible solution is to run
`seqkit` with the `rmdup` command to remove duplicates.

