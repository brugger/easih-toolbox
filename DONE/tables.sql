
CREATE TABLE runs (
  runfolder           VARCHAR(200) NOT NULL,
  platform            VARCHAR(100),
  stamp               BIGINT,
  status              VARCHAR(100) NOT NULL,

  PRIMARY KEY (runfolder,stamp)

);

CREATE TABLE offloads (
  runfolder           VARCHAR(200) NOT NULL,
  in_file              TEXT NOT NULL,
  out_file             TEXT NOT NULL,

  PRIMARY KEY (runfolder,out_file),
  KEY in_idx  (in_file(500)),
  KEY out_idx (out_file(500))
);


  
