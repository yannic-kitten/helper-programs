var error_info*: tuple[usage: string] = ("",)

proc usage_quit*(channel = stdout; error_code = QuitSuccess) =
  channel.writeLine("Usage: " & error_info.usage)
  error_code.quit()

proc error*(msg: string) =
  stderr.writeLine("Error: " & msg)
  stderr.usage_quit(QuitFailure)

proc error_if*(cond: bool; msg: string;) =
  if cond: error(msg)



# usage example:
when isMainModule:
  import std/cmdline
  error_info.usage = "./app arg1 arg2 arg3 ..."

  proc parseCli(args: seq[string]) =
    echo "Parsing CLI ..."
    # ...
    if "-h" in args:
      usage_quit()
    # ...
    if args.len < 3:
      error("missing arguments")
    error_if(args.len > 3, "too many arguments")
    # ...
    echo "CLI is valid"

  commandLineParams().parseCli
