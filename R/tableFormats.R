historyTableFormat <- c(
    'time               TEXT NOT NULL',
    'operation          TEXT NOT NULL',
    'MSsary_version     TEXT NOT NULL',
    'call               TEXT DEFAULT NULL',
    'augPackage         TEXT DEFAULT NULL',
    'augPackVersion     TEXT DEFAULT NULL'
)

# MsData

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
    'peakID     INTEGER PRIMARY KEY AUTOINCREMENT',
    'msLevel    INTEGER NOT NULL',
    'length     INTEGER NOT NULL',
    'mzMean     REAL NOT NULL',
    'maxHeight  REAL NOT NULL',
    'area       REAL NOT NULL',
    'peak       BLOB NOT NULL'
)

peakRtree <- c(
    'peakID',
    'scanStart',
    'scanEnd',
    'mzMin',
    'mzMax'
)

linkTableFormat <- c(
    'linkID     INTEGER PRIMARY KEY AUTOINCREMENT',
    'linkType   TEXT NOT NULL'
)

linkMemberFormat <- c(
    'linkID             INTEGER REFERENCES linkTable',
    'peakID             INTEGER REFERENCES peakInfo',
    'linkDescription    TEXT DEFAULT \'\''
)

# MsDataSet

memberTableFormat <- c(
    'memberID           INTEGER PRIMARY KEY AUTOINCREMENT',
    'saryLocation       TEXT NOT NULL',
    'memberName         TEXT NOT NULL',
    'referenceSample    INTEGER NOT NULL CHECK(referenceSample IN (0, 1)) DEFAULT 0'
)

groupTableFormat <- c(
    'groupID    INTEGER PRIMARY KEY AUTOINCREMENT',
    'msLevel    INTEGER NOT NULL',
    'nPeaks     INTEGER NOT NULL',
    'nMembers   INTEGER NOT NULL'
)

groupRtree <- c(
    'groupID',
    'rtStart',
    'rtEnd',
    'mzMin',
    'mzMax'
)

groupMemberFormat <- c(
    'groupID    INTEGER REFERENCES groupTable',
    'memberID   INTEGER REFERENCES members',
    'peakID     INTEGER NOT NULL'
)

rtCorTableFormat <- c(
    'memberID       INTEGER REFERENCES members',
    'scanNum        INTEGER NOT NULL',
    'retentionTime  REAL NOT NULL'
)