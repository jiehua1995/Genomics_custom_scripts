/*
 *Created Date: Tuesday, February 18th 2025, 12:36:43 pm
 *Author: Jie Hua
 *
 *Copyright (c) 2025 Jie Hua<Jie.Hua@lmu.de>
*/


#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define MAX_SEQ_LENGTH 10000

// Function to calculate GC content of a sequence
double calculate_gc_content(const char* sequence) {
    int gc_count = 0;
    int seq_len = strlen(sequence);
    
    for (int i = 0; i < seq_len; i++) {
        if (sequence[i] == 'G' || sequence[i] == 'C') {
            gc_count++;
        }
    }

    return (seq_len > 0) ? (double)gc_count / seq_len * 100 : 0;
}

// Function to calculate Nx statistic
int calculate_nx(int* lengths, int num_seqs, int x) {
    int threshold = 0;
    for (int i = 0; i < num_seqs; i++) {
        threshold += lengths[i];
    }
    threshold = threshold * x / 100;

    int cumulative_length = 0;
    for (int i = 0; i < num_seqs; i++) {
        cumulative_length += lengths[i];
        if (cumulative_length >= threshold) {
            return lengths[i];
        }
    }
    return 0;
}

// Function to process the FASTA file and print results
void process_fasta(const char* fasta_file, int machine_readable, int summarize, const char* output_file) {
    FILE* f = fopen(fasta_file, "r");
    if (!f) {
        fprintf(stderr, "Error opening file: %s\n", fasta_file);
        exit(1);
    }

    char line[1000];
    char seq_id[100];
    char sequence[MAX_SEQ_LENGTH];
    int seq_lengths[1000];
    double gc_contents[1000];
    int seq_count = 0;

    while (fgets(line, sizeof(line), f)) {
        // If line starts with '>', it's a sequence ID
        if (line[0] == '>') {
            if (seq_count > 0) {
                gc_contents[seq_count - 1] = calculate_gc_content(sequence);
                seq_lengths[seq_count - 1] = strlen(sequence);
            }
            // Save the sequence ID
            sscanf(line, ">%s", seq_id);
            memset(sequence, 0, sizeof(sequence));  // Reset sequence
        } else {
            // Add sequence data
            strcat(sequence, line);
        }
    }

    // Handle last sequence
    if (seq_count > 0) {
        gc_contents[seq_count - 1] = calculate_gc_content(sequence);
        seq_lengths[seq_count - 1] = strlen(sequence);
    }

    fclose(f);

    // Output results
    if (summarize) {
        // Sort the lengths in descending order
        for (int i = 0; i < seq_count - 1; i++) {
            for (int j = i + 1; j < seq_count; j++) {
                if (seq_lengths[i] < seq_lengths[j]) {
                    int temp = seq_lengths[i];
                    seq_lengths[i] = seq_lengths[j];
                    seq_lengths[j] = temp;
                }
            }
        }

        // Calculate N10, N20, N30, N50, N70, N80, N90
        printf("Summary Statistics:\n");
        printf("N10: %d\n", calculate_nx(seq_lengths, seq_count, 10));
        printf("N20: %d\n", calculate_nx(seq_lengths, seq_count, 20));
        printf("N30: %d\n", calculate_nx(seq_lengths, seq_count, 30));
        printf("N50: %d\n", calculate_nx(seq_lengths, seq_count, 50));
        printf("N70: %d\n", calculate_nx(seq_lengths, seq_count, 70));
        printf("N80: %d\n", calculate_nx(seq_lengths, seq_count, 80));
        printf("N90: %d\n", calculate_nx(seq_lengths, seq_count, 90));

        printf("Top 10 sequence lengths: ");
        for (int i = 0; i < 10 && i < seq_count; i++) {
            printf("%d ", seq_lengths[i]);
        }
        printf("\n");
    } else {
        // Print in human-readable or machine-readable format
        if (machine_readable) {
            printf("ID\tLength\tGC_Content\n");
            for (int i = 0; i < seq_count; i++) {
                printf("%s\t%d\t%.2f\n", seq_id, seq_lengths[i], gc_contents[i]);
            }
        } else {
            printf("%-20s%-10s%-15s\n", "ID", "Length", "GC_Content (%)");
            for (int i = 0; i < seq_count; i++) {
                printf("%-20s%-10d%-15.2f\n", seq_id, seq_lengths[i], gc_contents[i]);
            }
        }
    }
}

// Main function to parse arguments and call processing functions
int main(int argc, char *argv[]) {
    if (argc < 2) {
        fprintf(stderr, "Usage: %s <fasta_file> [-T] [-summarize] [-o <output_file>]\n", argv[0]);
        return 1;
    }

    int machine_readable = 0;
    int summarize = 0;
    const char* output_file = NULL;
    const char* fasta_file = argv[1];

    for (int i = 2; i < argc; i++) {
        if (strcmp(argv[i], "-T") == 0) {
            machine_readable = 1;
        } else if (strcmp(argv[i], "-summarize") == 0) {
            summarize = 1;
        } else if (strcmp(argv[i], "-o") == 0 && i + 1 < argc) {
            output_file = argv[i + 1];
            i++;
        }
    }

    process_fasta(fasta_file, machine_readable, summarize, output_file);
    return 0;
}
