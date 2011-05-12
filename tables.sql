#DROP TABLE snp;
#DROP TABLE reference;
#DROP TABLE population;
#DROP TABLE populations;
#DROP TABLE flags;

CREATE TABLE meta (
  version               VARCHAR(10),
  species               VARCHAR(120),
  reference		VARCHAR(120),  
  filedata              VARCHAR(30)

);


CREATE TABLE snp (

  rs                     VARCHAR(50),
  chr                    VARCHAR(6),
  pos			 INT NOT NULL,
  ref_base		 VARCHAR(100),
  alt_base		 VARCHAR(100),
  class                  VARCHAR(30),
  hgmd                   VARCHAR(1),
  flags			 VARCHAR(200),

  PRIMARY KEY (chr,pos),
  KEY rs_idx (rs)
);

CREATE TABLE flags (
  short VARCHAR(25) PRIMARY KEY NOT NULL,
  descr  VARCHAR(200)
);
