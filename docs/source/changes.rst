.. _changes:

Changes in Ringtail
######################
#TODO
2.0.0
Have to explicitly save receptor to calculate interactions (vina)

some cmdline keywords have changed (interactions mostly)

the logger now writes to file in addition to STDOUT

duplicate handling now works as expected, and is not a trait of the table. Must use the keyword in the call dealing with suspected duplicate results. If IGNORE it will not add results to results or inteeractions to bitvectors and interaction table. If REPLACE it will update the fields in results that change (not the unique columns), and update the bitvector fields, and for Interactions it will delete existing and insert new ones, maintaining pose_id but getting new table id.
A bitvector table and a table of indices for unique interactions are written once the database writing has been completed. These two tables will be remade every time the database is written to. This ensures duplicate handling constraints is handled accurately. this takes about 100 ms extra per 500 results.

An extra table has been added (Interactions) and a table has changed names (Interaction_bitvectors -> Interaction_bitvector_strings).

Method to update the database created, BUT If you created the database with the duplicate handling option, there is a chance of inconsistent behavior of anything involving interactions as
        the Pose_ID was not used as an explicit foreign key in db v1.0.0 and v1.1.0.
