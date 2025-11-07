import os
import hashlib
import argparse
from tqdm import tqdm  # Import tqdm for progress bar

def calculate_hash(file_path, hash_algorithm="sha256"):
    """
    Calculate the hash value of a given file.
    
    Parameters:
    - file_path: Path to the file for which we want to calculate the hash.
    - hash_algorithm: Hashing algorithm to use. Default is "sha256". 
      Other options include "md5", "sha1", etc.
    
    Returns:
    - A string representing the hash value of the file (hexadecimal format).
    """
    # Initialize the hash function using the specified algorithm (default: sha256)
    hash_func = hashlib.new(hash_algorithm)
    
    # Open the file in binary read mode
    with open(file_path, "rb") as f:
        # Read the file in chunks (8KB at a time) to handle large files
        while chunk := f.read(8192):  # Read 8KB at a time
            hash_func.update(chunk)  # Update the hash function with the chunk of data
    
    # Return the hexadecimal digest of the file's hash
    return hash_func.hexdigest()

def get_files_in_directory(directory):
    """
    Retrieve all file paths in the given directory (recursively).
    
    Parameters:
    - directory: The directory whose files we want to scan.
    
    Returns:
    - A list of file paths within the directory.
    """
    file_paths = []  # List to store the paths of all files
    
    # Walk through the directory tree, including subdirectories
    for root, _, files in os.walk(directory):
        # Add the full path of each file to the file_paths list
        for file in files:
            file_paths.append(os.path.join(root, file))
    
    return file_paths

def write_hashes_to_file(directory, output_file, hash_algorithm):
    """
    Calculate the hash of all files in the directory and write the file names and hashes to a text file.
    
    Parameters:
    - directory: The directory containing the files.
    - output_file: The name of the output text file where the results will be saved.
    - hash_algorithm: The hash algorithm to be used (e.g., 'sha256', 'md5', 'sha1').
    """
    # Get all file paths from the directory
    file_paths = get_files_in_directory(directory)
    
    # Open the output file in write mode
    with open(output_file, "w") as f:
        # Add a header to indicate the script source
        f.write("# This script is from: https://github.com/jiehua1995/Genomics_custom_scripts\n")
        f.write("# File Name\tHash Value\n")  # Optional header for clarity
        
        # Create a progress bar with tqdm
        for file_path in tqdm(file_paths, desc="Processing files", unit="file"):
            # Extract the file name from the file path (just the name, not the full path)
            file_name = os.path.basename(file_path)
            
            # Calculate the hash of the file using the specified algorithm
            file_hash = calculate_hash(file_path, hash_algorithm)
            
            # Write the file name and its corresponding hash to the output file
            f.write(f"{file_name}\t{file_hash}\n")

def main():
    # Set up argument parser
    parser = argparse.ArgumentParser(description="Calculate file hashes and output to a text file.")
    parser.add_argument("-p", "--path", required=True, help="Path to the directory to scan for files.")
    parser.add_argument("-o", "--output", required=True, help="Output file to save the hashes and file names.")
    parser.add_argument("-a", "--algorithm", default="sha256", choices=["md5", "sha1", "sha256", "sha224", "sha512"],
                        help="The hash algorithm to use. Default is 'sha256'. Options: 'md5', 'sha1', 'sha256', 'sha224', 'sha512'.")
    
    # Parse the command-line arguments
    args = parser.parse_args()

    # Call the function to process the directory and save results
    write_hashes_to_file(args.path, args.output, args.algorithm)

    # Notify the user that the hash values have been saved to the file
    print(f"Hash values have been saved to {args.output} using the {args.algorithm} algorithm.")

# Run the script if it is executed as a standalone program
if __name__ == "__main__":
    main()
