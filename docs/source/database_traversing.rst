.. _database_traversing:

Exploring the database in the command line
############################################
View the data contained within the database using a terminal, we recommend using the [VisiData tool](https://www.visidata.org/). In addition to command line visualization, this tool has a number of other feature, including ploting. Upon opening the database with `vd`, the terminal should look like this:

![Screenshot from 2022-05-18 14-57-22](https://user-images.githubusercontent.com/41704502/169162632-3a71d338-faa1-4109-8f04-40a96ee6d24e.png)

In this example (made with DLGs), the database contains ~3 poses for 9999 discrete ligands. Each of the rows here is a separate table or view within the database. From this screen, you can easily perform the sanity checks outline below. One should note that the number of column displayed on the first screen is 1 greater than the actual number of columns in a table (the number is correct for views). To more fully explore a given table, one may use the arrow keys or mouse to navigate to it, then press `Enter/Return` to access that table/view. The user may then scroll horizontally with the arrow keys, or press `q` to return up a level.

Using `vd` is particularly helpful to examine possible interactions of interest, stored within the `Interaction_indices` table.

To exit, return to the screen shown in the image above by pressing `q`, then press `q` to exit.

### Data integrity sanity checks
There are a few quick checks the user can make to ensure that the data has been properly written from the input files to the database. Discrepancies may indicate an error occurred while writing the database or the input file format did not match that which Ringtail expected.
- The number of rows in the `Ligands` table should match the number of input ligand files
- The number of rows in the `Results` and `Interaction_bitvectors` tables should match
- Number of columns in the `Interactions_bitvectors` table should match the number of rows in the `Interaction_indices` table + 1 (+2 if using `vd`)
- The number of rows in the `Results` table should be ~`max_poses`\* `number of files` and should be less than or equal to that number. For DLGs not every ligand may have up to `max_poses`, which is why the number of rows is typically smaller than `max_poses`\* `number of DLGs`.
- No ligand should have more than `max_poses` rows in the `Results` table.
- If storing all poses, the number of rows in the Results table should match the `number of ligands` * `number of output poses`.
