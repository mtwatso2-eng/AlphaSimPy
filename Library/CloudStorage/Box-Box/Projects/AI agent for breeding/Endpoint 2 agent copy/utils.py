import os
def list_folder_contents(folder_path):
        all_paths = []
        for root, dirs, files in os.walk(folder_path):
            for d in dirs:
                all_paths.append(os.path.join(root, d))
            for f in files:
                all_paths.append(os.path.join(root, f))
        return all_paths