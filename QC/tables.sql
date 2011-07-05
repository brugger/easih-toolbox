#DROP TABLE project;
#DROP TABLE sample;
#DROP TABLE file;
#DROP TABLE qv_boxplot;
#DROP TABLE qv_histogram;
#DROP TABLE base_distribution;
#DROP TABLE GC_distribution;
#DROP TABLE duplicates;
#DROP TABLE duplicated_seqs;
#DROP TABLE adaptors;
#DROP TABLE illumina_multiplex_stats;
#DROP TABLE illumina_lane_stats;
#DROP TABLE run;


CREATE TABLE project (

  pid                 INT NOT NULL AUTO_INCREMENT PRIMARY KEY,
  name                VARCHAR(3) NOT NULL ,

  KEY name_idx (name)

);



CREATE TABLE sample (

  sid                 INT NOT NULL AUTO_INCREMENT PRIMARY KEY,
  pid                 INT NOT NULL,
  name                VARCHAR(8) NOT NULL ,

  KEY name_idx (name),
  KEY pid_idx  (pid)
);


CREATE TABLE offloading (
  rid                 INT NOT NULL,
  stamp               BIGINT,
  status              VARCHAR(100) NOT NULL,

  PRIMARY KEY (rid,stamp)

);

CREATE TABLE illumina_lane_stats (

  rid                 INT NOT NULL,
  fid                 INT ,
  lane		      INT NOT NULL,
  read_nr	      INT NOT NULL,
  sample              VARCHAR(20) NOT NULL ,
  total_reads	      INT,
  pass_filter	      INT,

  PRIMARY KEY (rid, lane, read_nr),
  KEY rid_idx (rid),
  KEY fid_idx (fid)

);

CREATE TABLE illumina_multiplex_stats (
  rid                 INT NOT NULL,
  fid                 INT NOT NULL,
  lane		      INT NOT NULL,
  sample	      VARCHAR(20),
  bcode		      VARCHAR(20),
  ratio		      float,
  total_reads	      INT,
  pass_filter	      INT,

  PRIMARY KEY (rid, fid, lane, bcode),
  KEY rid_idx (rid),
  KEY fid_idx (fid)
  
);

CREATE TABLE run (

  rid                 INT NOT NULL AUTO_INCREMENT PRIMARY KEY,
  name                VARCHAR(100) NOT NULL ,
  platform            VARCHAR(50),

  KEY name_idx (name)
);


CREATE TABLE file (

  fid                 INT NOT NULL AUTO_INCREMENT PRIMARY KEY,
  sid                 INT NOT NULL ,
  name                VARCHAR(50) NOT NULL ,
  timestamp           BIGINT,
  rid		      INT NOT NULL,

  total_reads	      int,
  sample_size	      int,
  Q30bases	      float,
  duplicates	      float,
  partial_adaptors    float,
  Avg_AC	      float,
  
  KEY name_idx (name),
  KEY sid_idx  (sid)
);

CREATE TABLE mapping_stats (

  fid1                INT NOT NULL PRIMARY KEY,
  fid2                INT,
  reference           VARCHAR(100) NOT NULL,
  unique_hits	      INT,
  non_unique_hits     INT,
  duplicates          INT,
  
  KEY f12_idx (fid1, fid2)
);


CREATE TABLE qv_boxplot (
  fid                 INT NOT NULL,
  x            	      INT NOT NULL,
  q0		      float,
  q1		      float,
  q2		      float,
  q3		      float,
  q4		      float,
  
  PRIMARY KEY (fid, x),
  KEY fid_idx( fid )
);


CREATE TABLE qv_histogram (
  fid                 INT NOT NULL,
  x            	      INT NOT NULL,
  height	      float,
  
  PRIMARY KEY (fid, x),
  KEY fid_idx( fid )
);


CREATE TABLE base_distribution (
  fid                 INT NOT NULL,
  x            	      INT NOT NULL,
  A		      float,
  C		      float,
  G		      float,
  T		      float,
  N		      float,
  
  PRIMARY KEY (fid, x),
  KEY fid_idx( fid )
);


CREATE TABLE GC_distribution (
  fid                 INT NOT NULL,
  x            	      INT NOT NULL,
  percent_gc	      float,
  
  PRIMARY KEY (fid, x),
  KEY fid_idx( fid )
);


CREATE TABLE duplicates (
  fid                 INT NOT NULL,
  x            	      INT NOT NULL,
  observations	      INT,

  PRIMARY KEY (fid, x),
  KEY fid_idx( fid )
);

CREATE TABLE duplicated_seqs (
  fid                 INT NOT NULL,
  sequence	      varchar(400),
  percentage	      float,
  source	      varchar(100),
  
  PRIMARY KEY (fid, sequence),
  KEY fid_idx( fid )
);


CREATE TABLE adaptors (
  fid                 INT NOT NULL,
  x            	      INT NOT NULL,
  percent	      float,
  
  PRIMARY KEY (fid, x),
  KEY fid_idx( fid )
);




  
