import os


def folder_creator(options):
    """
    Creates a Folder in the selected directory

    :paran options: Name of the folder
    :type options: string

    returns: folder in the current directory
    """
    options = str(options.file_name)
    if not os.path.exists(options):
        os.mkdir(options)
        print("\nDirectory", options, "created at:\n",
              os.path.abspath(options))
    else:
        print("\nDirectory", options, "already exists")
