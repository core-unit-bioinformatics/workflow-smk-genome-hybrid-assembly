
localrules: copy_verkko_executables
rule copy_verkko_executables:
    output:
        exe = DIR_PROC.joinpath("tooling/verkko.exec.bak"),
        root_folder = directory(
            DIR_PROC.joinpath("tooling/verkko")
        ),
        lib_folder = directory(
            DIR_PROC.joinpath("tooling/verkko/lib"),
        )
    container:
        config["verkko_container"]
    shell:
        "mkdir -p {output.root_folder}"
            " && "
        "cp /repo/verkko/bin/verkko {output.exe}"
            " && "
        "mkdir -p {output.lib_folder}"
            " && "
        "cp -r /repo/verkko/lib/ {output.root_folder}"
            " && "
        "chmod u+x {output.exe}"
            " && "
        "chmod -R u+x {output.root_folder}/lib/verkko/bin"


localrules: copy_verkko_testdata
rule copy_verkko_testdata:
    output:
        hifi = DIR_PROC.joinpath("testdata/verkko/hifi.fastq.gz"),
        nano = DIR_PROC.joinpath("testdata/verkko/ont.fastq.gz"),
    container:
        config["verkko_container"]
    shell:
        "mkdir -p proc/testdata/verkko"
            " && "
        "cp /testdata/hifi.fastq.gz {output.hifi}"
            " && "
        "cp /testdata/ont.fastq.gz {output.nano}"


localrules: adapt_verkko_launcher
rule adapt_verkko_launcher:
    input:
        source = DIR_PROC.joinpath("tooling/verkko.exec.bak"),
    output:
        target = DIR_PROC.joinpath("tooling/verkko/bin/verkko"),
    run:
        import io
        import os
        import stat
        inject_lines = {
            "add_variable": "snake_hpc_profile=\"\"\n",
            "add_hpc_switch": ""
                + "    "
                + "elif [ \"$opt\" = \"--hpc\" ] ;"
                + "                then grid=\"hpc\";\n",
            "add_hpc_profile": ""
                + "    "
                + "elif [ \"$opt\" = \"--hpc-profile\" ] ;"
                + "        then snake_hpc_profile=$arg; shift\n",
            "add_hpc_help": ""
                + "    "
                + "echo \"    --hpc"
                + "                    "
                + "Enable generic HPC/scheduler (no built-in support).\"\n",
            "add_profile_help": ""
                + "    "
                + "echo \"    --hpc-profile"
                + "            "
                + "For generic HPC/scheduler use, set path to Snakemake profile.\"\n",
            "set_hpc_profile": ""
                + "elif [ $grid = \"hpc\" ] ; then\n"
                + "    "
                + "echo >> ${outd}/snakemake.sh \""
                + "  --profile ${snake_hpc_profile} \\\\\"\n"
                + "    "
                + "echo >> ${outd}/snakemake.sh \"  --local-cores 1 \\\\\"\n"
        }

        new_script = io.StringIO()
        profile_trigger = False
        with open(input.source, "r") as script:
            for ln, line in enumerate(script):
                if profile_trigger and line.startswith("else"):
                    add_section = inject_lines["set_hpc_profile"]
                    new_script.write(add_section)
                    new_script.write(line)
                    continue
                else:
                    if "if [ $grid = \"local\" ]" in line:
                        profile_trigger = True
                    new_script.write(line)
                if line.startswith("snakeopts="):
                    add_line = inject_lines["add_variable"]
                    new_script.write(add_line)
                if "[ \"$opt\" = \"--lsf\" ]" in line:
                    add_line = inject_lines["add_hpc_switch"]
                    new_script.write(add_line)
                    add_line = inject_lines["add_hpc_profile"]
                    new_script.write(add_line)
                if "Enable IBM Spectrum LSF support." in line:
                    add_line = inject_lines["add_hpc_help"]
                    new_script.write(add_line)
                    add_line = inject_lines["add_profile_help"]
                    new_script.write(add_line)

        with open(output.target, "w") as script:
            _ = script.write(new_script.getvalue())
        os.chmod(
            output.target,
            stat.S_IXUSR | stat.S_IWRITE | stat.S_IREAD | stat.S_IRGRP
        )
    # END OF RUN BLOCK


rule run_verkko_test_local:
    input:
        exe = rules.adapt_verkko_launcher.output.target,
        hifi = rules.copy_verkko_testdata.output.hifi,
        nano = rules.copy_verkko_testdata.output.nano,
    output:
        wd = directory(
            DIR_PROC.joinpath("testdata/verkko/local/assembly")
        )
    log:
        DIR_LOG.joinpath("testdata/verkko/local/assembly.log")
    benchmark:
        DIR_RSRC.joinpath("testdata/verkko/local/assembly.log")
    conda:
        "../../envs/verkko_env.yaml"
    threads: 4
    resources:
        mem_mb = lambda wildcards, attempt: 12288 * attempt,
        mem_gb = lambda wildcards, attempt: 12 * attempt,
        time_hrs = lambda wildcards, attempt: attempt * attempt,
    params:
        bin_folder = lambda wc, input: str(input.exe).strip("/verkko")
    shell:
        "export PATH=$PATH:$PWD/{params.bin_folder} ; "
        "/usr/bin/time -v "
        "verkko --local "
            "--local-memory {resources.mem_gb} "
            "--local-cpus {threads} "
            "--threads {threads} "
            "--python `which python` "
            "--mbg `which MBG` "
            "--graphaligner `which GraphAligner` "
            "--hifi {input.hifi} "
            "--nano {input.nano} "
            "-d {output.wd} &> {log}"


localrules: run_verkko_test_cluster
rule run_verkko_test_cluster:
    input:
        exe = rules.adapt_verkko_launcher.output.target,
        hifi = rules.copy_verkko_testdata.output.hifi,
        nano = rules.copy_verkko_testdata.output.nano,
        profile = config["verkko_smk_profile"],
    output:
        wd = directory(
            DIR_PROC.joinpath("testdata/verkko/cluster/assembly")
        )
    log:
        DIR_LOG.joinpath("testdata/verkko/cluster/assembly.log")
    benchmark:
        DIR_RSRC.joinpath("testdata/verkko/cluster/assembly.log")
    conda:
        "../../envs/verkko_env.yaml"
    params:
        bin_folder = lambda wc, input: str(input.exe).strip("/verkko")
    shell:
        "export PATH=$PATH:$PWD/{params.bin_folder} ; "
        "/usr/bin/time -v "
        "verkko --hpc --hpc-profile $PWD/{input.profile} "
            "--python `which python` "
            "--mbg `which MBG` "
            "--graphaligner `which GraphAligner` "
            "--hifi {input.hifi} "
            "--nano {input.nano} "
            "-d {output.wd} &> {log}"


localrules: run_verkko_tests
rule run_verkko_tests:
    input:
        assm_dirs = [
            DIR_PROC.joinpath("testdata/verkko/local/assembly"),
            DIR_PROC.joinpath("testdata/verkko/cluster/assembly")
        ]
