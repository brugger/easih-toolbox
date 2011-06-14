
CREATE TABLE runs (
  runfolder           VARCHAR(200) NOT NULL,
  platform            VARCHAR(100),
  ts                  TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
  status              VARCHAR(100) NOT NULL,

  PRIMARY KEY (runfolder,ts,status)

);
