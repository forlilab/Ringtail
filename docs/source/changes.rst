.. _changes:

Changes in Ringtail
######################
#TODO
2.0.0
Have to explicitly save receptor to calculate interactions (vina)

some cmdline keywords have changed (interactions mostly)

the logger now writes to file in addition to STDOUT

duplicate handling now works as expected, and is not a trait of the table. Must use the keyword in the call dealing with suspected duplicate results. If IGNORE it will not add results to results or inteeractions to bitvectors and interaction table. If REPLACE it will update the fields in results that change (not the unique columns), and update the bitvector fields, and for Interactions it will delete existing and insert new ones, maintaining pose_id but getting new table id.
bitvector table no longer created, the information in that table is created in-situ for the two filter methods that use it, using pandas