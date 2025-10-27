import tkinter as tk
from tkinter import filedialog, messagebox, ttk
import os
import threading
from collections import Counter

class FastaAnalyzer:
    def __init__(self):
        self.root = tk.Tk()
        self.root.title("FASTA File Analyzer")
        self.root.geometry("800x600")
        self.create_widgets()
        
    def create_widgets(self):
        main_frame = ttk.Frame(self.root, padding="10")
        main_frame.grid(row=0, column=0, sticky=(tk.W, tk.E, tk.N, tk.S))
        
        self.root.columnconfigure(0, weight=1)
        self.root.rowconfigure(0, weight=1)
        main_frame.columnconfigure(1, weight=1)
        main_frame.rowconfigure(2, weight=1)
        
        ttk.Label(main_frame, text="Select FASTA File:").grid(row=0, column=0, sticky=tk.W, pady=5)
        self.file_path_var = tk.StringVar()
        self.file_entry = ttk.Entry(main_frame, textvariable=self.file_path_var, state="readonly", width=50)
        self.file_entry.grid(row=0, column=1, sticky=(tk.W, tk.E), padx=5, pady=5)
        
        self.browse_button = ttk.Button(main_frame, text="Browse", command=self.browse_file)
        self.browse_button.grid(row=0, column=2, padx=5, pady=5)
        
        self.analyze_button = ttk.Button(main_frame, text="Analyze FASTA File", command=self.analyze_file, state="disabled")
        self.analyze_button.grid(row=1, column=0, columnspan=3, pady=10)
        
        self.progress_var = tk.StringVar(value="Ready")
        ttk.Label(main_frame, textvariable=self.progress_var).grid(row=2, column=0, columnspan=3, sticky=tk.W, pady=5)
        
        self.progress_bar = ttk.Progressbar(main_frame, mode='indeterminate')
        self.progress_bar.grid(row=3, column=0, columnspan=3, sticky=(tk.W, tk.E), pady=5)
        
        ttk.Label(main_frame, text="Analysis Results:").grid(row=4, column=0, sticky=tk.W, pady=(10, 5))
        
        text_frame = ttk.Frame(main_frame)
        text_frame.grid(row=5, column=0, columnspan=3, sticky=(tk.W, tk.E, tk.N, tk.S), pady=5)
        text_frame.columnconfigure(0, weight=1)
        text_frame.rowconfigure(0, weight=1)
        
        self.results_text = tk.Text(text_frame, wrap=tk.WORD, height=20, width=80)
        scrollbar = ttk.Scrollbar(text_frame, orient=tk.VERTICAL, command=self.results_text.yview)
        self.results_text.configure(yscrollcommand=scrollbar.set)
        
        self.results_text.grid(row=0, column=0, sticky=(tk.W, tk.E, tk.N, tk.S))
        scrollbar.grid(row=0, column=1, sticky=(tk.N, tk.S))
        
    def browse_file(self):
        file_path = filedialog.askopenfilename(
            title="Select FASTA File",
            filetypes=[
                ("FASTA files", "*.fasta *.fa *.fna *.ffn *.frn *.faa"),
                ("All files", "*.*")
            ]
        )
        
        if file_path:
            self.file_path_var.set(file_path)
            self.analyze_button.config(state="normal")
            
            try:
                file_size = os.path.getsize(file_path)
                size_mb = file_size / (1024 * 1024)
                self.progress_var.set(f"Selected: {os.path.basename(file_path)} ({size_mb:.2f} MB)")
            except OSError:
                self.progress_var.set(f"Selected: {os.path.basename(file_path)}")
    
    def analyze_file(self):
        file_path = self.file_path_var.get()
        if not file_path:
            messagebox.showerror("Error", "Please select a FASTA file first.")
            return
        
        self.analyze_button.config(state="disabled")
        self.progress_bar.start()
        self.progress_var.set("Analyzing FASTA file...")
        
        thread = threading.Thread(target=self._analyze_fasta_thread, args=(file_path,))
        thread.daemon = True
        thread.start()
    
    def _analyze_fasta_thread(self, file_path):
        try:
            results = self.parse_fasta_file(file_path)
            self.root.after(0, self._display_results, results)
        except Exception as e:
            self.root.after(0, self._show_error, str(e))
    
    def parse_fasta_file(self, file_path):
        results = {
            'file_path': file_path,
            'file_size': os.path.getsize(file_path),
            'sequences': [],
            'total_sequences': 0,
            'total_length': 0,
            'overall_composition': Counter()
        }
        
        current_header = None
        current_sequence = []
        
        with open(file_path, 'r', encoding='utf-8', errors='replace') as file:
            for line in file:
                line = line.strip()
                
                if line.startswith('>'):
                    if current_header and current_sequence:
                        seq_data = self._process_sequence(current_header, ''.join(current_sequence))
                        results['sequences'].append(seq_data)
                        results['total_length'] += seq_data['length']
                        results['overall_composition'].update(seq_data['composition'])
                    
                    current_header = line[1:]
                    current_sequence = []
                else:
                    if current_header:
                        current_sequence.append(line)
        
        if current_header and current_sequence:
            seq_data = self._process_sequence(current_header, ''.join(current_sequence))
            results['sequences'].append(seq_data)
            results['total_length'] += seq_data['length']
            results['overall_composition'].update(seq_data['composition'])
        
        results['total_sequences'] = len(results['sequences'])
        return results
    
    def _process_sequence(self, header, sequence):
        composition = Counter(sequence.upper())
        length = len(sequence)
        
        lines = sequence.split('\n')
        lines_80_chars = all(len(line) == 80 for line in lines[:-1]) and len(lines[-1]) <= 80
        
        return {
            'header': header,
            'sequence': sequence,
            'length': length,
            'composition': composition,
            'lines_80_chars': lines_80_chars,
            'gc_content': self._calculate_gc_content(composition)
        }
    
    def _calculate_gc_content(self, composition):
        total = sum(composition.values())
        if total == 0:
            return 0.0
        gc_count = composition.get('G', 0) + composition.get('C', 0)
        return (gc_count / total) * 100
    
    def _display_results(self, results):
        self.progress_bar.stop()
        self.analyze_button.config(state="normal")
        self.progress_var.set("Analysis complete")
        
        self.results_text.delete(1.0, tk.END)
        
        output = []
        output.append("=" * 80)
        output.append("FASTA FILE ANALYSIS RESULTS")
        output.append("=" * 80)
        output.append(f"File: {os.path.basename(results['file_path'])}")
        output.append(f"File size: {results['file_size'] / (1024*1024):.2f} MB")
        output.append(f"Number of sequences: {results['total_sequences']}")
        output.append(f"Total sequence length: {results['total_length']:,} bases")
        output.append("")
        
        if results['total_length'] > 0:
            output.append("OVERALL COMPOSITION:")
            output.append("-" * 40)
            for base in sorted(results['overall_composition'].keys()):
                count = results['overall_composition'][base]
                percentage = (count / results['total_length']) * 100
                output.append(f"{base}: {count:,} ({percentage:.2f}%)")
            
            gc_count = results['overall_composition'].get('G', 0) + results['overall_composition'].get('C', 0)
            gc_percentage = (gc_count / results['total_length']) * 100
            output.append(f"GC content: {gc_count:,} ({gc_percentage:.2f}%)")
            output.append("")
        
        for i, seq in enumerate(results['sequences'], 1):
            output.append(f"SEQUENCE {i}:")
            output.append(f"  Header: {seq['header']}")
            output.append(f"  Length: {seq['length']:,} bases")
            output.append(f"  Lines 80 chars: {'Yes' if seq['lines_80_chars'] else 'No'}")
            output.append(f"  GC content: {seq['gc_content']:.2f}%")
            output.append("  Composition:")
            for base in sorted(seq['composition'].keys()):
                count = seq['composition'][base]
                percentage = (count / seq['length']) * 100
                output.append(f"    {base}: {count:,} ({percentage:.2f}%)")
            output.append("")
        
        self.results_text.insert(tk.END, '\n'.join(output))
        self.results_text.see(tk.END)
    
    def _show_error(self, error_message):
        self.progress_bar.stop()
        self.analyze_button.config(state="normal")
        self.progress_var.set("Error occurred")
        messagebox.showerror("Analysis Error", f"An error occurred during analysis:\n{error_message}")
    
    def run(self):
        self.root.mainloop()

def main():
    app = FastaAnalyzer()
    app.run()

if __name__ == "__main__":
    main()