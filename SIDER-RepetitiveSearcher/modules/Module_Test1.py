import os

def folder_creator(Folder_Name):
    """
    Creates a Folder in the selected directory
    """
    if not os.path.exists(Folder_Name):
        os.mkdir(Folder_Name)
        print("\nDirectory", Folder_Name, "created at:\n",
             os.path.abspath(Folder_Name))
    else:
        print("\nDirectory", Folder_Name, "already exists")