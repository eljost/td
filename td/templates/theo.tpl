    <!doctype html>
<html>
    <head>
        <meta charset="utf-8">
        <style>
            html * {
                font-size: 1.1em;
            }
            .state {
                margin-bottom: 2%;
            }
            .mo {
                height: 15%;
                min-height: 10em;
                width: 15%;
                min-width: 10em;
                border: 1px solid #000000;
                padding: 1%;
            }
            sub {
                font-size: .8em;
            }
        </style>
    </head>
    <body>
    {% for state in states %}
    <div class="state">
        <h2>
            S<sub>{{ state.id }}</sub>,
            {{ state.spat }},
            {{ "%.1f" | format(state.l) }} nm,
            {{ "%.2f" | format(state.dE) }} eV,
            f={{ "%.4f" | format(state.f) }}
        </h2>
        {% for i in range((nto_dict[state.id]|length // 2)+1) %}
        <div>
            <figure>
                <img class="mo" src="theodore/{{ nto_dict[state.id][i][0] }}" />
                <img class="mo" src="theodore/{{ nto_dict[state.id][i][1] }}" />
                <figcaption>{{ "%.1f" | format(nto_dict[state.id][i][2]*100) }}%</figcaption>
            </figure>
        </div>
        {% endfor %}
        </div>
        {% endfor %}
    </body>
</html>
