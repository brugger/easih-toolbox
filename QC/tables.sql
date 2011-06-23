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


CREATE TABLE file (

  fid                 INT NOT NULL AUTO_INCREMENT PRIMARY KEY,
  sid                 INT NOT NULL ,
  name                VARCHAR(50) NOT NULL ,
  platform            VARCHAR(50),
  timestamp           BIGINT,

  sample_size	      int,
  Q30bases	      float,
  duplicates	      float,
  partial_adaptors    float,
  Avg_AC	      float,
  


  KEY name_idx (name),
  KEY sid_idx  (sid)
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
  
  PRIMARY KEY (fid, x),
  KEY fid_idx( fid )
);


CREATE TABLE adaptors (
  fid                 INT NOT NULL,
  x            	      INT NOT NULL,
  percent	      float,
  
  PRIMARY KEY (fid, x),
  KEY fid_idx( fid )
);




  
