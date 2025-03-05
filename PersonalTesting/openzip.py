import zipfile
import os
zipFilePath = "sdf.zip"
extractDir = "sdf_files"
if os.path.exists(zipFilePath):

    with zipfile.ZipFile(zipFilePath, 'r') as zipRef:
        zipRef.extractall(extractDir)
        print(f"Files extracted to: {extractDir}")
else:
    print(f"ZIP file {zipFilePath} not found.")
