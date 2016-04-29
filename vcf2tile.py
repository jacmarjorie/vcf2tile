import sys, os, vcf, uuid, json, subprocess, tarfile

from collections import OrderedDict


CONST_TILEDB_FIELDS = OrderedDict()
CONST_TILEDB_FIELDS["END"]             = { "vcf_field_class" : ["INFO"],          "type": "int" }
CONST_TILEDB_FIELDS["BaseQRankSum"]    = { "vcf_field_class" : ["INFO"],          "type":"float" }
CONST_TILEDB_FIELDS["ClippingRankSum"] = { "vcf_field_class" : ["INFO"],          "type":"float" }
CONST_TILEDB_FIELDS["MQRankSum"]       = { "vcf_field_class" : ["INFO"],          "type":"float" }
CONST_TILEDB_FIELDS["ReadPosRankSum"]  = { "vcf_field_class" : ["INFO"],          "type":"float" }
CONST_TILEDB_FIELDS["MQ"]              = { "vcf_field_class" : ["INFO"],          "type":"float" }
CONST_TILEDB_FIELDS["MQ0"]             = { "vcf_field_class" : ["INFO"],          "type":"int"   }
CONST_TILEDB_FIELDS["AF"]              = { "vcf_field_class" : ["INFO"],          "type":"float", "length":"A" }
CONST_TILEDB_FIELDS["AN"]              = { "vcf_field_class" : ["INFO"],          "type":"int",   "length":1 }
CONST_TILEDB_FIELDS["AC"]              = { "vcf_field_class" : ["INFO"],          "type":"int",   "length":"A" }
CONST_TILEDB_FIELDS["DP"]              = { "vcf_field_class" : ["INFO","FORMAT"], "type":"int"   }
CONST_TILEDB_FIELDS["MIN_DP"]          = { "vcf_field_class" : ["FORMAT"],        "type":"int"   }
CONST_TILEDB_FIELDS["GQ"]              = { "vcf_field_class" : ["FORMAT"],        "type":"int"   }
CONST_TILEDB_FIELDS["SB"]              = { "vcf_field_class" : ["FORMAT"],        "type":"int",   "length":4 }
CONST_TILEDB_FIELDS["AD"]              = { "vcf_field_class" : ["FORMAT"],        "type":"int",   "length":"R" }
CONST_TILEDB_FIELDS["PL"]              = { "vcf_field_class" : ["FORMAT"],        "type":"int",   "length":"G" }
CONST_TILEDB_FIELDS["GT"]              = { "vcf_field_class" : ["FORMAT"],        "type":"int",   "length":"P" }
CONST_TILEDB_FIELDS["PS"]              = { "vcf_field_class" : ["FORMAT"],        "type":"int",   "length":1 }


def writeVIDMappingFile(reader, vid_map_file, fields_dict=CONST_TILEDB_FIELDS):
    vid_mapping = OrderedDict()
    vid_mapping["fields"] = fields_dict
    vid_mapping["contigs"] = OrderedDict()
    contigs = vid_mapping["contigs"]
    offset = 0
    # use vcf header to construct contigs
    for reference in reader.contigs:
      contigs[reference] = {"length": reader.contigs[reference].length,
                                 "tiledb_column_offset": offset}
      offset += (long(reader.contigs[reference].length) + 1000)

    writeJSON2File(vid_mapping, vid_map_file)


def getCallSets(reader, filename, callset_check, sampleTag = False, row_counter=0):

  callsets = OrderedDict()

  for i in range(0, len(reader.samples)):
    if sampleTag is False:
      name = reader.samples[i]
      idx_in_file = i
    else:
      name = reader.metadata['SAMPLE'][i]['SampleName']
      idx_in_file = reader.samples.index(reader.metadata['SAMPLE'][i]['ID'])

    # check duplicate
    if name in callset_check:
      nuuid = str(uuid.uuid4())
      sys.stderr.write('Duplicate callset name '+name+' : appending _'+nuuid+'\n');
      name += ('_'+nuuid)

    callsets[name] = OrderedDict()
    callsets[name]['row_idx'] = row_counter
    callsets[name]['idx_in_file'] = idx_in_file
    callsets[name]['filename'] = filename 
    print callsets[name]
    print 'callset', name, 'with', str(row_counter)

    row_counter += 1


  return callsets, row_counter


def writeJSON2File(input_json, output_file):
  with open(output_file, "w") as outFP:
    json.dump(input_json, outFP, indent=2, separators=(',', ': '))


if __name__ == "__main__":
  import argparse 
  parser = argparse.ArgumentParser(description = "Load a list of VCFs into TileDB.")

  # includes path to vid_map and callset_map
  # workspace, array
  parser.add_argument("-l", "--loader", required=True, type=str,
                      help="base loader file")
  parser.add_argument("-i", "--input", type=str, required=True, help="VCF file to be imported.")
  parser.add_argument("-s", "--sampleTag", action="store_true", required=False, help="Use SAMPLE tag in VCF header to name callsets.")
  parser.add_argument("-L", "--load", action="store_true", required=False, help="Run vcf2tiledb after creating necessary configs.")
  parser.add_argument("-t", "--tar", action="store_true", required=False, help="Extract tar.")  
  args = parser.parse_args()


  # check if the file is a tarball
  # this is specific to the VC use case
  # data prep
  inputs = []
  if args.tar:
    tar = tarfile.open(args.input, 'r')
    for member in tar.getmembers():
      if member.name[-6:] == 'vcf.gz':
        inputs.append(tar.extractfile(member.name))

  else:
    inputs.append(open(args.input, 'rb'))
  
  with open(args.loader) as conf:
    config = json.load(conf)

  callset_map_file = config['callset_mapping_file']
  vid_map_file = config['vid_mapping_file']

  # if there are callsets in the array
  if os.path.isfile(callset_map_file):
    nuuid = str(uuid.uuid4())
    with open(callset_map_file) as callset_file:
      callset_mapping = json.load(callset_file)
      callsets = callset_mapping['callsets']
      # set new loading position
      last = callsets[max(callsets, key=lambda v:callsets[v]['row_idx'])]['row_idx'] + 1
  # else, make the file and start fresh
  else:
    callset_mapping = {}
    callset_mapping['callsets'] = {}
    callsets = callset_mapping['callsets']


  # update callsets, write vid if new
  rc = config.get('lb_callset_row_idx', 0)
  for input_file in inputs:

    r = vcf.Reader(input_file)
    # set schema from the first vcf that is imported
    # if this file exists, then it won't be wrote
    if not os.path.isfile(vid_map_file):
      writeVIDMappingFile(r, vid_map_file)

    new_callset, rc = getCallSets(r, input_file.name, callsets, sampleTag=args.sampleTag, row_counter=rc)
    callsets.update(new_callset)

      
  writeJSON2File(callset_mapping, callset_map_file)
  config['lb_callset_row_idx'] = rc
  writeJSON2File(config, args.loader)
	
  if args.load:
  	processArgs = ['vcf2tiledb', os.path.abspath(args.loader)]
 	# load the files
  	pipe = subprocess.Popen(processArgs, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
  	output, error = pipe.communicate()

  	if pipe.returncode != 0:
    		raise Exception("subprocess run: {0}\nFailed with stdout: \n-- \n{1} \n--\nstderr: \n--\n{2} \n--".format(" ".join(processArgs), output, error))

  for inp in inputs:
    inp.close()
