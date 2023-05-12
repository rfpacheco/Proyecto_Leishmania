First it reads a .csv file and takes the rows which include in the second column (the one with the chromosome IDs), the defined chromosome.
Then, it creates a .csv file  with all the selectes rows in a subfolder in a specified directory using :func:`~modules.files_manager.folder_creator` and :func:`~modules.files_manager.csv_creator`.

.. admonition:: Example of use

   If ``chromosome_ID = LinJ.02``, it will take all rows in the **.csv** where the second column contains ``LinJ.02``

:param path_input: Path to the main CSV file where data will be filtered. Initially is a CSV file wich was output from :func:`~modules.blaster.blastn_blaster` alone to obtain the initial data.
:type path_input: string

:param chromosome_ID: Chromosome identifier, e.g., *LinJ.07*, which were present more then once in the ``path_input``. See :func:`~modules.blaster.repetitive_blaster`for more information.
:type chromosome_ID: member of a list

:param main_folder_path: Path where the results will be placed. It will create a subfolder with the ``chromosome_ID`` name.
:type main_folder_path: string

:return: A .csv file with the selected "chromosome_ID".
:rtype: csv file
