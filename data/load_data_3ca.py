from dataclasses import dataclass
from pathlib import Path
import pandas as pd
from scipy.io import mmread
import numpy as np
from typing import Optional, Union
from scipy import sparse


@dataclass
class DATA3CA:
    folder_path: Union[str, Path]
    counts_file: Union[str, Path] = 'Exp_data_UMIcounts.mtx'
    cell_data_file: Union[str, Path] = 'Cells.csv'
    gene_data_file: Union[str, Path] = 'Genes.txt'
    meta_data_file: Union[str, Path] = 'Meta-data.csv'

    _counts: Optional[np.ndarray] = None
    _cells: Optional[pd.DataFrame] = None
    _genes: Optional[list] = None
    _metadata: Optional[pd.DataFrame] = None

    def __post_init__(self):
        self.folder_path = Path(self.folder_path)
        self.counts_file = self.folder_path / self.counts_file
        self.cell_data_file = self.folder_path / self.cell_data_file
        self.gene_data_file = self.folder_path / self.gene_data_file
        self.meta_data_file = self.folder_path / self.meta_data_file

        self._validate_files()

    def _validate_files(self) -> None:
        """Validate existence of all required files."""
        for file_path in [self.counts_file, self.cell_data_file,
                          self.gene_data_file, self.meta_data_file]:
            if not file_path.exists():
                raise FileNotFoundError(f"File {file_path} does not exist")

    @property
    def counts(self) -> np.ndarray:
        """Lazy load counts data."""
        if self._counts is None:
            self._counts = mmread(self.counts_file).transpose()
        if not sparse.isspmatrix_csr(self._counts):
            self._counts = self._counts.tocsr()
        return self._counts

    @property
    def cells(self) -> pd.DataFrame:
        """Lazy load cell data."""
        if self._cells is None:
            self._cells = pd.read_csv(self.cell_data_file, header=0, index_col=0)
        return self._cells

    @property
    def genes(self) -> list:
        """Lazy load gene data."""
        if self._genes is None:
            with open(self.gene_data_file) as f:
                self._genes = [line.strip() for line in f]
        return self._genes

    @property
    def metadata(self) -> pd.DataFrame:
        """Lazy load metadata."""
        if self._metadata is None:
            self._metadata = pd.read_csv(self.meta_data_file, header=0, index_col=0)
        return self._metadata

    def validate_data_consistency(self) -> None:
        """Validate data dimensions match across files."""
        counts_shape = self.counts.shape
        assert counts_shape[0] == len(self.cells), "Cell count mismatch"
        assert counts_shape[1] == len(self.genes), "Gene count mismatch"
