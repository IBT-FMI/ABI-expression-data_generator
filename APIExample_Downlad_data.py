import argparse
import copy
import json
import os
import sys
import shutil
import urllib
import urllib.request
import zipfile
import numpy
import nibabel
import re
import socket
import tarfile
from socket import timeout
from collections import defaultdict
from mhd_utils_3d import *
from nipype.interfaces.ants import ApplyTransforms
#from nipype.interfaces.fsl import fslorient

API_SERVER = "http://api.brain-map.org/"
API_DATA_PATH = API_SERVER + "api/v2/data/"

#GLOBAL DEFAULT = None, so no timeout error. Try this??
socket.setdefaulttimeout(100)

def GetGeneNames(startRow=0,numRows=2000,totalRows=-1):
    """
    Queries the Allen Mouse Brain Institute website for all gene expression data available for download.

    Returns:
    --------
    info: defaultdict
        key: genename, value: list of corresponding SectionDataSetID (SectionDataSet: see "http://help.brain-map.org/display/api/Data+Model")
        ID needed to specify download target.

    """

    startRow = startRow
    numRows = numRows
    totalRows = totalRows
    rows = []
    GeneNames = []
    SectionDataSetID = []
    info = defaultdict(list)

    done = False

    while not done:
        pagedUrl = API_DATA_PATH +"query.json?criteria=model::SectionDataSet,rma::criteria,[failed$eqfalse],products[abbreviation$eq'Mouse'],treatments[name$eq'ISH'],rma::include,genes,specimen(donor(age)),plane_of_section" + '&startRow=%d&numRows=%d' % (startRow,numRows)

        print(pagedUrl)
        source = urllib.request.urlopen(pagedUrl).read()
        response = json.loads(source)
        rows += response['msg']
        for x in response['msg']:

            if x['failed'] == False:
            #if x['failed'] == False and x['expression'] == True :
                info[x['genes'][0]['acronym']].append(x['id'])

        if totalRows < 0:
            totalRows = int(response['total_rows'])

        startRow += len(response['msg'])

        if startRow >= totalRows:
            done = True

    return info

def download_all_ISH(info,name="ABI_geneexpression_data",version="9999"):
    """
    Downloads all datasets corresponding to SectionDataSetID given, converts data format from mhd/raw to nii and transforms data to dsurqec template.

    Parameters:
   -----------
    info: defaultdict
        key: genename, value: list of corresponding SectionDataSetID (SectionDataSet: see "http://help.brain-map.org/display/api/Data+Model")
        ID needed to specify download target.

    """
    failed_downloads = list()
    data_path = name + "-" + version
    #TODO: script keeps hanging somewhere, maybe timeout for connection, or try catch block and saving exp files numbers for later downloads:
    if not os.path.exists(data_path): os.mkdir(data_path)
    download_url = "http://api.brain-map.org/grid_data/download/"
    for gene in info.items():
        gene_name = gene[0]
        gene_ids = gene[1]
        #replace brackets with '_' and remove ',*
        #TODO:put into sep. fcn to avoid doing it twice
        gene_r = re.sub('[()]',"_",gene_name)
        #gene_r = re.sub('\W', '',gene_r)  dont do that, will replace '-' as well
        gene_r = re.sub("'","",gene_r)
        gene_r = re.sub('\*',"",gene_r)
        path_to_gene = os.path.join("/mnt/data/setinadata/abi_data/geneexpression/ABI_geneexpression_data-9999",gene_r)
        #TODO: change logic, check if right file exists
        if os.path.exists(path_to_gene):continue
        if not os.path.exists(path_to_gene) : os.mkdir(path_to_gene)
        print(gene_r)
        for id in gene_ids:
            url = download_url + str(id)
            try:
                fh = urllib.request.urlretrieve(url)
            except timeout:
                print("timeout with " + str(id)+ gene_r)
                failed_downloads.append(id)
                shutil.rmtree(path_to_gene)
                continue

            zf = zipfile.ZipFile(fh[0])
            filename = str.split((fh[1]._headers[6][1]),'filename=')[1]
            filename = str.split(filename,'.zip')[0]
            #replace brackets with '_' and remove ',*
            filename = re.sub('[()]',"_",filename)
            filename = re.sub('\*', '',filename)
            filename = re.sub("'","",filename)
            path_to_folder = os.path.join(path_to_gene,filename)
            zf.extractall(os.path.join(path_to_gene,filename))
            zf.close()

            #some datasets without energy.mhd file. Skip and delete folder
            if not os.path.isfile(os.path.join(path_to_folder,"energy.mhd")):
                    print("removing" + str(id)+gene_r)
                    shutil.rmtree(path_to_folder)
                    continue

            path_to_mhd = os.path.join(path_to_folder,"energy.mhd")
            path_to_nifti = convert_raw_to_nii(path_to_mhd,filename)
            apply_composite(path_to_nifti)
            os.remove(path_to_nifti)

    if len(failed_downloads) > 0:
        print("failed: ")
        for item in failed_downloads:
            print(str(item))

def convert_raw_to_nii(input_file,output_file):
    """
    Converts mhd/raw format to NIfTI and orients data matrix in RAS-space.

    Parameters:
    -----------
        input_file : str
            path to .mhd file
        output_file : str
            filename prefix

    Returns:
    ---------
        output_path : str
            path to generated NIfTI - file
    """
    path = os.path.abspath('.')
    image_array, meta_header = load_raw_data_with_mhd(input_file)

    #Read header infomormation and create affine matrix
    dims = numpy.array(meta_header["ElementSpacing"].split(" "),dtype=numpy.float)
    affine_matrix = numpy.zeros((4,4),dtype=numpy.float)
    affine_matrix[0,0] = dims[0]
    affine_matrix[1,1] = dims[1]
    affine_matrix[2,2] = dims[2]

    #Orient in RAS
    image_array = numpy.swapaxes(image_array,1,2)

    image_array = image_array[::-1,:,:]
    image_array = image_array[:,:,::-1]
    image_array = image_array[:,::-1,:]

    #Bring to the right units
    affine_matrix = affine_matrix*0.001

    affine_matrix[3,3] = 1

    img = nibabel.Nifti1Image(image_array,affine_matrix)
    name = output_file + '.nii.gz'
    output_path = os.path.join(os.path.dirname(input_file),name)
    nibabel.save(img,output_path)
    
    return output_path

def apply_composite(file):
    """
    Uses ANTS ApplyTransforms to transform image to 

    Parameters :
    ------------

    file : str
        path to image

    """
    at = ApplyTransforms()
    at.inputs.dimension = 3
    at.inputs.input_image = file
    at.inputs.reference_image = '/usr/share/mouse-brain-atlases/dsurqec_200micron_masked.nii'
    name = str.split(os.path.basename(file),'.nii')[0] + '_2dsurqec.nii.gz'
    at.inputs.interpolation = 'NearestNeighbor' #TODO: Sure??
    at.inputs.output_image = os.path.join(os.path.dirname(file),name)
    at.inputs.transforms = '/usr/share/mouse-brain-atlases/abi2dsurqec_Composite.h5'
    at.run()

    #TODO sform to qform

def save_info(info):
    f = open("ABI_genes_ids.csv","w")
    for gene in info:
        f.write('\n')
        f.write(gene)
        for id in info[gene]:
            f.write("," + str(id))

def create_archive(name,version):
   path = "ABI_geneexpression_data-9999"
   tar_name = name + "-" + version + ".tar.xz"
   with tarfile.open(tar_name, "w:xz") as tar_handle:
      for root,dirs,files in os.walk(path):
         for file in files:
            print(file)
            tar_handle.add(os.path.join(root,file))


def main():
   parser = argparse.ArgumentParser(description="Similarity",formatter_class=argparse.ArgumentDefaultsHelpFormatter)
   parser.add_argument('--package_name','-n',type=str,default="ABI_geneexpression_data")
   parser.add_argument('--package_version','-v',type=str,default="9999")
   parser.add_argument('--startRow','-s',type=int,default=0)
   parser.add_argument('--numRows','-r',type=int,default=2000)
   parser.add_argument('--totalRows','-t',type=int,default=-1)
   args=parser.parse_args()

   #info=GetGeneNames(startRow=args.startRow,numRows=args.numRows,totalRows=args.totalRows)
   #download_all_ISH(info,name=args.package_name,version=args.package_version)
   #save_info(info)
   create_archive(args.package_name,args.package_version)

if __name__ == "__main__":
    main()
