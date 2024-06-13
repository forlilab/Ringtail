.. _changes:

Changes in Ringtail
######################
#TODO
2.0.0
Have to explicitly save receptor to calculate interactions (vina)

some cmdline keywords have changed (interactions mostly)

the logger now writes to file in addition to STDOUT

Duplicate handling now works as expected, and is not a trait of the table. Must use the keyword in the call dealing with suspected duplicate results. If IGNORE it will not add results to results or inteeractions to bitvectors and interaction table. If REPLACE it will update the fields in results that change (not the unique columns), and update the bitvector fields, and for Interactions it will delete existing and insert new ones, maintaining pose_id but getting new table id.
An Interaction_bitvectors table and a table of Interaction_indices for unique interactions are written once the database writing has been completed. These two tables will be remade every time the database is written to. This ensures duplicate handling constraints is handled accurately. this takes about 100 ms extra per 500 results, or less than 5 minutes extra time for 1,000,000 files. 

An extra table has been added (Interactions) which allows for proper handling if user chooses to deal with duplicated entries in a specific way. This table also ensures 1-to-1 of the primary key in the Results table to a foreign key dealing with interactions. 

Method to update the database created, BUT If you created the database with the duplicate handling option, there is a chance of inconsistent behavior of anything involving interactions as the Pose_ID was not used as an explicit foreign key in db v1.0.0 and v1.1.0.
