import tkinter as tk
from tkinter import filedialog, messagebox, ttk
import os
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from collections import Counter

class SlidingWindowApp:
    def __init__(self):
        self.root = tk.Tk()
        self.root.title("FASTA Sliding Window Analyzer")
        self.root.geometry("1000x700")
        
        self.sequence = ""
        self.file_path = ""
        
        self.create_widgets()
        
    def create_widgets(self):
        main_frame = ttk.Frame(self.root, padding="10")
        main_frame.pack(fill=tk.BOTH, expand=True)
        
        ttk.Label(main_frame, text="Select FASTA File:").grid(row=0, column=0, sticky=tk.W, pady=5)
        self.file_path_var = tk.StringVar(value="No file selected")
        self.file_entry = ttk.Entry(main_frame, textvariable=self.file_path_var, state="readonly", width=50)
        self.file_entry.grid(row=0, column=1, sticky=(tk.W, tk.E), padx=5, pady=5)
        
        self.browse_button = ttk.Button(main_frame, text="Browse", command=self.browse_file)
        self.browse_button.grid(row=0, column=2, padx=5, pady=5)
        
        ttk.Label(main_frame, text="Window Size:").grid(row=0, column=3, sticky=tk.W, padx=(20,5), pady=5)
        self.window_var = tk.IntVar(value=30)
        self.window_entry = ttk.Spinbox(main_frame, from_=1, to=1000, textvariable=self.window_var, width=8)
        self.window_entry.grid(row=0, column=4, padx=5, pady=5)
        
        self.analyze_button = ttk.Button(main_frame, text="Analyze & Plot", command=self.analyze_sequence)
        self.analyze_button.grid(row=0, column=5, padx=15, pady=5)
        
        self.info_var = tk.StringVar(value="Load a FASTA file to begin.")
        ttk.Label(main_frame, textvariable=self.info_var).grid(row=1, column=0, columnspan=6, sticky=tk.W, pady=5)
        
        self.fig, self.ax = plt.subplots(figsize=(10, 6))
        self.canvas = FigureCanvasTkAgg(self.fig, main_frame)
        self.canvas.get_tk_widget().grid(row=2, column=0, columnspan=6, sticky=(tk.W, tk.E, tk.N, tk.S), pady=10)
        
        main_frame.columnconfigure(1, weight=1)
        main_frame.rowconfigure(2, weight=1)
        
    def browse_file(self):
        file_path = filedialog.askopenfilename(
            title="Select FASTA File",
            filetypes=[("FASTA files", "*.fasta *.fa *.fna *.ffn *.frn *.faa"), ("All files", "*.*")]
        )
        
        if file_path:
            try:
                with open(file_path, 'r', encoding='utf-8') as file:
                    sequence = ""
                    for line in file:
                        line = line.strip()
                        if not line.startswith('>'):
                            sequence += line.upper()
                
                if not sequence:
                    raise ValueError("No sequence data found")
                
                self.sequence = sequence
                self.file_path = file_path
                self.file_path_var.set(os.path.basename(file_path))
                self.info_var.set(f"Loaded sequence: {len(sequence)} bases")
                
            except Exception as e:
                messagebox.showerror("Error", f"Failed to load file: {str(e)}")
    
    def analyze_sequence(self):
        if not self.sequence:
            messagebox.showwarning("No sequence", "Please load a FASTA file first.")
            return
        
        try:
            window_size = self.window_var.get()
            if window_size <= 0:
                raise ValueError("Window size must be positive")
            
            positions, frequencies = self.sliding_window_analysis(self.sequence, window_size)
            
            if not positions:
                messagebox.showwarning("Invalid", f"Sequence length ({len(self.sequence)}) is smaller than window size ({window_size})")
                return
            
            self.plot_results(positions, frequencies)
            self.info_var.set(f"Analyzed {len(positions)} windows of size {window_size}")
            
        except Exception as e:
            messagebox.showerror("Error", f"Analysis failed: {str(e)}")
    
    def sliding_window_analysis(self, sequence, window_size):
        positions = []
        a_freqs = []
        c_freqs = []
        g_freqs = []
        t_freqs = []
        
        for i in range(len(sequence) - window_size + 1):
            window = sequence[i:i + window_size]
            total = len(window)
            
            if total == 0:
                continue
            
            a_count = window.count('A')
            c_count = window.count('C')
            g_count = window.count('G')
            t_count = window.count('T')
            
            positions.append(i + 1)
            a_freqs.append(a_count / total)
            c_freqs.append(c_count / total)
            g_freqs.append(g_count / total)
            t_freqs.append(t_count / total)
        
        return positions, {'A': a_freqs, 'C': c_freqs, 'G': g_freqs, 'T': t_freqs}
    
    def plot_results(self, positions, frequencies):
        self.ax.clear()
        
        self.ax.plot(positions, frequencies['A'], label='A', linewidth=2)
        self.ax.plot(positions, frequencies['C'], label='C', linewidth=2)
        self.ax.plot(positions, frequencies['G'], label='G', linewidth=2)
        self.ax.plot(positions, frequencies['T'], label='T', linewidth=2)
        
        self.ax.set_title("Sliding Window Nucleotide Frequencies")
        self.ax.set_xlabel("Window Position")
        self.ax.set_ylabel("Relative Frequency")
        self.ax.set_ylim(0, 1)
        self.ax.legend()
        self.ax.grid(True, alpha=0.3)
        
        self.canvas.draw()
    
    def run(self):
        self.root.mainloop()

def main():
    app = SlidingWindowApp()
    app.run()

if __name__ == "__main__":
    main()
