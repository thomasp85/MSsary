headerTableFormat <- c(
    'seqNum                     INTEGER NOT NULL',
    'acquisitionNum             INTEGER PRIMARY KEY',
    'msLevel                    INTEGER NOT NULL',
    'peaksCount                 INTEGER NOT NULL',
    'totIonCurrent              REAL NOT NULL',
    'retentionTime              REAL NOT NULL',
    'basePeakMZ                 REAL NOT NULL',
    'basePeakIntensity          REAL NOT NULL',
    'collisionEnergy            REAL',
    'ionisationEnergy           REAL',
    'lowMZ                      REAL NOT NULL',
    'highMZ                     REAL NOT NULL',
    'precursorScanNum           INTEGER',
    'precursorMZ                REAL',
    'precursorCharge            INTEGER',
    'precursorIntensity         REAL',
    'mergedScan                 INTEGER',
    'mergedResultScanNum        INTEGER',
    'mergedResultStartScanNum   INTEGER',
    'mergedResultEndScanNum     INTEGER'
)

scanTableFormat <- c(
    'scanNum            INTEGER UNIQUE FOREIGN KEY(scanNum) REFERENCES header(acquisitionNum)',
    'peaksCount         INTEGER NOT NULL',
    'totIonCurrent      REAL NOT NULL',
    'basePeakMZ         REAL NOT NULL',
    'basePeakIntensity  REAL NOT NULL',
    'lowMZ              REAL NOT NULL',
    'highMZ             REAL NOT NULL',
    'remove             INTEGER NOT NULL CHECK (remove IN (0, 1)) DEFAULT 0',
    'scan               BLOB'
)

peakTableFormat <- c(
    'peakID     INTEGER PRIMARY KEY',
    'scanStart  INTEGER FOREIGN KEY(scanStart) REFERENCES header(acquisitionNum)',
    'scanEnd    INTEGER FOREIGN KEY(scanEnd) REFERENCES header(acquisitionNum)',
    'mzMin      REAL NOT NULL',
    'mzMean     REAL NOT NULL',
    'mzMax      REAL NOT NULL',
    'maxHeight  REAL NOT NULL',
    'area       REAL NOT NULL',
    'peak       BLOB NOT NULL'
)