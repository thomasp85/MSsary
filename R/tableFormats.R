historyTableFormat <- c(
    'time               TEXT NOT NULL',
    'operation          TEXT NOT NULL',
    'MSsary_version     TEXT NOT NULL',
    'call               TEXT DEFAULT NULL',
    'augPackage         TEXT DEFAULT NULL',
    'augPackVersion     TEXT DEFAULT NULL'
)

headerTableFormat <- c(
    'seqNum                     INTEGER NOT NULL UNIQUE',
    'acquisitionNum             INTEGER PRIMARY KEY',
    'msLevel                    INTEGER NOT NULL',
    'polarity                   INTEGER NOT NULL',
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
    'scanNum            INTEGER UNIQUE REFERENCES header',
    'peaksCount         INTEGER',
    'totIonCurrent      REAL',
    'basePeakMZ         REAL',
    'basePeakIntensity  REAL',
    'lowMZ              REAL',
    'highMZ             REAL',
    'remove             INTEGER NOT NULL CHECK (remove IN (0, 1)) DEFAULT 0',
    'scan               BLOB'
)

peakTableFormat <- c(
    'peakID     INTEGER PRIMARY KEY',
    'scanStart  INTEGER REFERENCES header(seqNum)',
    'scanEnd    INTEGER REFERENCES header(seqNum)',
    'length     INTEGER NOT NULL',
    'mzMin      REAL NOT NULL',
    'mzMean     REAL NOT NULL',
    'mzMax      REAL NOT NULL',
    'maxHeight  REAL NOT NULL',
    'area       REAL NOT NULL',
    'peak       BLOB NOT NULL'
)