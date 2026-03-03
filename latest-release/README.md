# Directory Structure
.
├── README.md
├── piv_build
│   ├── bin
│   ├── lib
│   ├── libexec
│   └── share
├── piv_pipeline.py
└── tests
|   ├── pipeline_utility_testing
|   └── SVSout_23227179_1724441851
│       ├── jpgframes  
|       └── pivframes
└── testbatches
    ├── tmp.2YNbHiPCwK.txt
    ├── tmp.6PlqXRWBtN.txt
    .
    .
    .

# Test Instructions
This workflow does need working versions of Julia 1 and Python3. All packages for the Julia environment are wrapped up into the `piv_build` directory, editing this directory in any way is not recommended. Packages used in `piv_pipeline`.py are all included in the base installation of Python3.
 
**1.** `cd` into the root directory of this file. 
**2.** Change the permissions on `piv_pipeline.py` using `chmod +x piv_pipeline.py` and run the file using `./piv_pipeline.py` **or** skip updating the permissions and simply run `python3 piv_pipeline.py`.
**3.** MacOS doesn't like files downloaded from the internet, so run `sudo xattr -r -d com.apple.quarantine .` to allow files in this directory to run.
**4.** Run `./piv_pipeline.py -h` to show CLI options for the program. *NOTE: `--out` and `--in_path` are required for the script to run.*
**5.** Adjust the number of processes you'd like to spawn by editing the `NPROCS` constant at the top of `piv_pipeline.py`.
**6.** To run the test: `./piv_pipeline.py --out tests/pipeline_utility_testing/SVSout_23227179_1724441851/pivframes --in_path tests/pipeline_utility_testing/testbatches/ `.
**7.** .mat files are output to `pivframes/`.