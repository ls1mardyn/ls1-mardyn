import argparse
import json
import numpy as np

from validation_createJSON import parse_resultwriter_file, parse_log_file

def compare_data(new_data, validation_data, reltolerance):
    '''
    Compares new and validation data sets using numpy's isclose with a specified tolerance.
    '''
    differences = []
    for i, (new_entry, validation_entry) in enumerate(zip(new_data, validation_data)):
        entry_differences = {}
        for key in new_entry.keys():
            if key in validation_entry:
                if not np.isclose(new_entry[key], validation_entry[key], rtol=reltolerance):
                    entry_differences[key] = {
                        'presentRun': new_entry[key],
                        'reference': validation_entry[key]
                    }
        if entry_differences:
            if 'Simstep' in new_entry.keys():  # Log file
                simstep = new_entry['Simstep']
            elif 'simstep' in new_entry.keys():  # ResultWriter
                simstep = new_entry['simstep']
            else:
                simstep = np.nan
            differences.append({'index': i, 'simstep': simstep, 'differences': entry_differences})
    return differences

def compare_validation_file(validation_file, new_log_file):
    '''
    Compares the new files with the data stored in the validation JSON file.
    '''
    try:
        with open(validation_file, 'r') as f:
            validation_data = json.load(f)
    except FileNotFoundError:
        print(f'Error: Validation file "{validation_file}" not found.')
        exit(1)
    except json.JSONDecodeError:
        print(f'Error: Validation file "{validation_file}" is not a valid JSON file.')
        exit(1)
    except Exception as e:
        print(f'Failed with exception: {e}')
        raise

    # Relative tolerance; chosen so that small deviations due to number of ranks are neglected
    # Specified in metadata
    reltolerance = validation_data['metadata']['reltolerance']

    # Parse files and compare data; errors are handled in respective function

    # Process log file
    try:
        validation_log_data = validation_data['logfile']
    except KeyError:
        print('Error: Validation file is missing required key ("logfile").')
        exit(1)
    new_log_data = parse_log_file(new_log_file)
    log_diffs = compare_data(new_log_data, validation_log_data, reltolerance)

    # Process file of ResultWriter
    try:
        validation_result_data = validation_data['ResultWriter']
    except KeyError:
        print('Error: Validation file is missing required keys ("ResultWriter").')
        exit(1)

    # Get filename of output of ResultWriter from metadata
    # Gives "None" if no ResultWriter file was specified during generation of validation file
    new_result_file = validation_data['metadata']['ResultWriter_filename']

    if new_result_file == "None":
        new_result_data = []
    else:
        new_result_data = parse_resultwriter_file(new_result_file)
    
    resultwriter_diffs = compare_data(new_result_data, validation_result_data, reltolerance)
    
    return {
        'log_diffs': log_diffs,
        'resultwriter_diffs': resultwriter_diffs
    }

if __name__ == '__main__':
    '''
    Compares the output (Logger, ResultWriter) of a simulation using numpy's isclose function.
    Since ResultWriter writes to file but Logger to stdout, the logfile has to be specified.
    '''
    parser = argparse.ArgumentParser(
        description='Compare new simulation with a JSON validation file',
        epilog='Example usage: validation_compare_files.py --validation-file=validation.json --logfile=new_log.log --resultfile=new_result.res',
    )
    parser.add_argument('--validation-file', required=True, help='Path to the JSON validation file')
    parser.add_argument('--logfile', required=True, help='Path to the new log file')

    args = parser.parse_args()

    differences = compare_validation_file(args.validation_file, args.logfile)

    # Print differences
    if differences['log_diffs'] or differences['resultwriter_diffs']:
        print('Differences found:')
        print(json.dumps(differences, indent=4))
        exit(1)  # Exit with failure
    else:
        print('No differences found :-)')
        exit(0)  # Exit with success

