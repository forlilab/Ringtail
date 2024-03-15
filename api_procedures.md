# API procedures

### Creating a new ringtail core

A ringtail core is created by instantiating the object with a database. Currently, a database can only be added upon instantiation.

```
rtc = RingtailCore("output.db")
```

Options regarding logging level and storage properties automatically sets to default values (see table), and can be changed using methods described below. 
*In writing this I recognize that it is currently impossible to debug the instantiation of the object, unless default logger level is set to debug. 

### Adding results
To add results to the database, use the `add_results_from_files` method that takes as input files and directories to upload,
as well as a receptor path and database properties and how to handle the resutls (how many poses to save, how to deal with interactions if vina).
Files can be added either directly in the method call (`file`, `file_list` etc), or you can create a file object first, and use `file_source_object`.
The same can be done for the storage options/properties, by providing an `optionsdict` instead of individual options. 

```
rtc.add_results_from_files( file = None, 
                            file_path = None, 
                            file_list = None, 
                            file_pattern = None, 
                            recursive = None, 
                            receptor_file=None,
                            save_receptor = None
                            store_all_poses = None,
                            max_poses = None,
                            add_interactions = None,
                            interaction_tolerance = None,
                            interaction_cutoffs = None,
                            max_proc = None)
```
or a file sources object can be created separately, as well as the result processing/write options:
```

file_sources = rtc.set_file_sources(file_path = "files/directory/",
                                file_pattern = "*.dlg*", 
                                recursive = True, 
                                receptor_file = "receptor.pdbqt",
                                save_receptor = True)

writeoptions = rtc.set_results_processing_options(store_all_poses = None,
                                                    max_poses = None,
                                                    add_interactions = None,
                                                    interaction_tolerance = None,
                                                    interaction_cutoffs = None,
                                                    max_proc = None)


rtc.add_results_from_files(file_source_object = file_sources,
                            optionsdict = writeoptions)
```

--> need to change how files are handled a little bit, don't like how it is now
--> future, add interactions after the fact as an API method? 

### Filtering on a database 

[describe how filters are set and used]
- might want to add filters to the filter method? 
- or better more verbose to make filter object? Should be same for all

### Options
All options from Ringtail are organized in five different objects:
- General options: Docking mode, whether or not summary should be printed when a database is inquired, and logging level
- Storage options: 
[list all options, values, method to set them]

- need to reorganize and document the options file, add typehints to core, and ensure description is same cmd line and ringtail options (and in ringtailcore hinting)

### Using the config file