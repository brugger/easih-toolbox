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
#DROP TABLE illumina_sample_stats;
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
  total_reads	      INT,
  pass_filter	      INT,
  total_bases	      INT,
  QV30_bases	      INT,

  PRIMARY KEY (rid, lane),
  KEY rid_idx (rid),

);

CREATE TABLE illumina_sample_stats (
  rid                 INT NOT NULL,
  fid                 INT NOT NULL,
  lane		      INT NOT NULL,
  read_nr	      INT NOT NULL,
  sample	      VARCHAR(20),
  bcode		      VARCHAR(20),
  ratio		      float,
  total_reads	      INT,

  PRIMARY KEY (rid, fid, lane, read_nr, bcode),
  KEY rid_idx (rid),
  KEY fid_idx (fid)
  
);


CREATE TABLE instrument (

  mid                 INT NOT NULL AUTO_INCREMENT PRIMARY KEY,
  name                VARCHAR(100) NOT NULL UNIQUE,

  KEY name_idx (name)
);

CREATE TABLE run (

  rid                 INT NOT NULL AUTO_INCREMENT PRIMARY KEY,
  name                VARCHAR(100) NOT NULL UNIQUE,
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
  read_length         int,
  sample_size	      int,
  Q30bases	      float,
  duplicates	      float,
  partial_adaptors    float,
  Avg_AC	      float,
  
  KEY name_idx (name),
  KEY sid_idx  (sid)
);

CREATE TABLE mapping_stats (

  fid1                INT NOT NULL,
  fid2                INT,
  reference           VARCHAR(100) NOT NULL,
  unique_hits	      INT,
  non_unique_hits     INT,
  duplicates          INT,
  
  PRIMARY KEY (fid1, reference),
  KEY f1_idx (fid1),
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
  q5		      float,
  q6		      float,
  q7		      float,
  
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


DELETE FROM project;
DELETE FROM sample;
DELETE FROM file;
DELETE FROM qv_boxplot;
DELETE FROM qv_histogram;
DELETE FROM base_distribution;
DELETE FROM GC_distribution;
DELETE FROM duplicates;
DELETE FROM duplicated_seqs;
DELETE FROM adaptors;
DELETE FROM illumina_multiplex_stats;
DELETE FROM illumina_lane_stats;
DELETE FROM run;
