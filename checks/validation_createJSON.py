import argparse
import json
import re
import os
from datetime import datetime
import pandas as pd

def parse_log_file(filepath):
    '''
    Parses the output of the Logger.

    Args:
        filepath (str): Path to the log file.

    Returns:
        list[dict]: A list of dictionaries representing the parsed data.
    '''
    try:
        data = []
        # '\d.eE+-' to also match scientific notation and positive/negative numbers
        pattern = re.compile(r'Simstep = (\d+)\s+T = ([\d.eE+-]+)\s+U_pot = ([\d.eE+-]+)\s+p = ([\d.eE+-]+)')

        with open(filepath, 'r') as f:
            for line in f:
                match = pattern.search(line)
                if match:
                    data.append({
                        'Simstep': int(match.group(1)),
                        'T': float(match.group(2)),
                        'U_pot': float(match.group(3)),
                        'p': float(match.group(4))
                    })
        return data
    except Exception as e:
        raise ValueError(f'Failed to parse log file: {filepath}. Error: {e}')

def parse_resultwriter_file(filepath):
    '''
    Parses the file from the ResultWriter using pandas to handle variable columns

    Args:
        filepath (str): Path to the result file.

    Returns:
        list[dict]: A list of dictionaries representing the parsed data.
    '''
    try:
        # Read the file, skipping comments and ignoring lines starting with '#'
        df = pd.read_csv(filepath, delim_whitespace=True, comment='#', engine='python')
        return df.to_dict(orient='records')
    except Exception as e:
        raise ValueError(f'Failed to parse file from ResultWriter: {filepath}. Error: {e}')

def create_validation_file(log_data, result_data, output_file):
    '''
    Combines the data into a single JSON file with metadata.

    Args:
        result_data (dict with name(string) and data(list[dict])): Parsed data from ResultWriter
        log_data (list[dict]): Parsed data from output of Logger
        output_file (str): Path to the output JSON file.
    '''
    
    try:
        commithash = os.popen('git rev-parse --short HEAD').read().strip()
    except:
        commithash = ''
    
    validation_data = {
        'metadata': {
            'created_at': datetime.now().isoformat(),
            'commit_at_creation': commithash,
            'build_options': 'CC=gcc CXX=g++ cmake -DVECTOR_INSTRUCTIONS=AVX2 -DCMAKE_BUILD_TYPE=Release -DENABLE_AUTOPAS=OFF -DAUTOPAS_ENABLE_RULES_BASED_AND_FUZZY_TUNING=ON -DENABLE_ALLLBL=OFF -DOPENMP=ON -DENABLE_MPI=ON -DENABLE_UNIT_TESTS=ON -DENABLE_VTK=ON ..',
            'comment': '',
            'reltolerance': 1e-8,  # Relative tolerance used for comparison
            'ResultWriter_filename': result_data['name'],
        },
        'logfile' : log_data,
        'ResultWriter' : result_data['data'],
    }

    with open(output_file, 'w') as f:
        json.dump(validation_data, f, indent=4)

    print(f'Validation file created: {output_file}')

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Process ls1 output files and create a JSON file used for validation',
        epilog='Example usage: validation_createJSON.py --logfile=out.log --resultfile=result.res --output=validation.json'
    )
    parser.add_argument('--logfile', required=True, help='Path to the log file')
    parser.add_argument('--resultfile', help='Path to the file created by the ResultWriter')
    parser.add_argument('--output', default='validation.json', help='Path to the output JSON validation file')

    args = parser.parse_args()

    result_data = dict()

    # Parse the input files
    if args.resultfile is not None:
        result_data['data'] = parse_resultwriter_file(args.resultfile)
        result_data['name'] = os.path.basename(args.resultfile)
    else:
        result_data['data'] = list(dict())
        result_data['name'] = 'None'
    log_data = parse_log_file(args.logfile)

    # Create the validation file
    create_validation_file(log_data, result_data, args.output)
