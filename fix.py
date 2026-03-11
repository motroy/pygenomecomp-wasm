import re

with open('index.html', 'r') as f:
    html = f.read()

with open('pygenomecomp_wasm.py', 'r') as f:
    py_code = f.read()

# Find the pymodule script tag block and replace its contents
start_tag = '<script type="text/python" id="pymodule">\n'
end_tag = '\n</script>'

start_idx = html.find(start_tag)
if start_idx != -1:
    end_idx = html.find(end_tag, start_idx)
    if end_idx != -1:
        new_html = html[:start_idx + len(start_tag)] + py_code + html[end_idx:]
        with open('index.html', 'w') as f:
            f.write(new_html)
        print("Updated index.html")
    else:
        print("Could not find end tag")
else:
    print("Could not find start tag")
